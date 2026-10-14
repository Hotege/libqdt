#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <future>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "qdt.h"
#include "wrapper.hxx"

namespace parallel
{
    class ThreadPool
    {
    public:
        explicit ThreadPool(size_t n) : stop_(false)
        {
            if (n == 0) n = 1;
            for (size_t i = 0; i < n; ++i)
                workers_.emplace_back([this] { worker_loop(); });
        }

        ~ThreadPool()
        {
            {
                std::unique_lock<std::mutex> lock(mtx_);
                stop_ = true;
            }
            cv_.notify_all();
            for (auto& t : workers_) t.join();
        }

        template <class F>
        auto enqueue(F&& f) -> std::future<decltype(f())>
        {
            using R = decltype(f());
            auto task = std::make_shared<std::packaged_task<R()>>(
                std::forward<F>(f));
            std::future<R> fut = task->get_future();
            {
                std::unique_lock<std::mutex> lock(mtx_);
                tasks_.emplace([task]() { (*task)(); });
            }
            cv_.notify_one();
            return fut;
        }

        bool try_run_one()
        {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(mtx_);
                if (tasks_.empty()) return false;
                task = std::move(tasks_.front());
                tasks_.pop();
            }
            task();
            return true;
        }

        template <class T>
        T wait_for(std::future<T>& fut)
        {
            for (;;)
            {
                if (fut.wait_for(std::chrono::seconds(0)) ==
                    std::future_status::ready)
                    return fut.get();
                if (!try_run_one())
                    return fut.get();
            }
        }

    private:
        void worker_loop()
        {
            for (;;)
            {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(mtx_);
                    cv_.wait(lock, [this] { return stop_ || !tasks_.empty(); });
                    if (stop_ && tasks_.empty()) return;
                    task = std::move(tasks_.front());
                    tasks_.pop();
                }
                task();
            }
        }

        std::vector<std::thread> workers_;
        std::queue<std::function<void()>> tasks_;
        std::mutex mtx_;
        std::condition_variable cv_;
        bool stop_;
    };

    std::atomic<int> g_next_thread_id{ 0 };
    thread_local int tls_thread_id = -1;

    inline int get_thread_id()
    {
        if (tls_thread_id < 0)
            tls_thread_id = g_next_thread_id.fetch_add(1, std::memory_order_relaxed);
        return tls_thread_id;
    }

    using Point = qdt_vec2d;

    namespace robust
    {
        namespace exact
        {
            struct Exp
            {
                static constexpr int LOCAL_CAP = 32;
                double local_buf[LOCAL_CAP];
                double* heap_ptr = nullptr;
                int size = 0;
                int capacity = LOCAL_CAP;

                Exp() = default;
                ~Exp() { if (heap_ptr) delete[] heap_ptr; }

                Exp(const Exp& other)
                {
                    size = other.size;
                    if (other.heap_ptr != nullptr)
                    {
                        capacity = other.capacity;
                        heap_ptr = new double[capacity];
                        std::copy(other.heap_ptr, other.heap_ptr + size, heap_ptr);
                    }
                    else
                    {
                        capacity = LOCAL_CAP;
                        std::copy(other.local_buf, other.local_buf + size, local_buf);
                    }
                }

                Exp(Exp&& other) noexcept
                {
                    size = other.size;
                    capacity = other.capacity;
                    heap_ptr = other.heap_ptr;
                    if (!heap_ptr)
                        std::copy(other.local_buf, other.local_buf + size, local_buf);
                    other.heap_ptr = nullptr;
                    other.size = 0;
                    other.capacity = LOCAL_CAP;
                }

                Exp& operator=(const Exp& other)
                {
                    if (this == &other) return *this;
                    if (other.heap_ptr != nullptr)
                    {
                        if (capacity < other.capacity)
                        {
                            if (heap_ptr) delete[] heap_ptr;
                            capacity = other.capacity;
                            heap_ptr = new double[capacity];
                        }
                        else if (!heap_ptr)
                        {
                            heap_ptr = new double[other.capacity];
                            capacity = other.capacity;
                        }
                        std::copy(other.heap_ptr, other.heap_ptr + other.size, heap_ptr);
                    }
                    else
                    {
                        if (heap_ptr) { delete[] heap_ptr; heap_ptr = nullptr; }
                        capacity = LOCAL_CAP;
                        std::copy(other.local_buf, other.local_buf + other.size, local_buf);
                    }
                    size = other.size;
                    return *this;
                }

                Exp& operator=(Exp&& other) noexcept
                {
                    if (this == &other) return *this;
                    if (heap_ptr) delete[] heap_ptr;
                    size = other.size;
                    capacity = other.capacity;
                    heap_ptr = other.heap_ptr;
                    if (!heap_ptr)
                        std::copy(other.local_buf, other.local_buf + size, local_buf);
                    other.heap_ptr = nullptr;
                    other.size = 0;
                    other.capacity = LOCAL_CAP;
                    return *this;
                }

                inline void clear() { size = 0; }
                inline double* data() { return heap_ptr ? heap_ptr : local_buf; }
                inline const double* data() const { return heap_ptr ? heap_ptr : local_buf; }

                inline void push_back(double val)
                {
                    if (size >= capacity)
                    {
                        int new_cap = capacity * 2;
                        double* new_heap = new double[new_cap];
                        std::copy(data(), data() + size, new_heap);
                        if (heap_ptr) delete[] heap_ptr;
                        heap_ptr = new_heap;
                        capacity = new_cap;
                    }
                    data()[size++] = val;
                }

                inline bool empty() const { return size == 0; }
                inline double back() const { return size > 0 ? data()[size - 1] : 0.0; }
                inline double operator[](int i) const { return data()[i]; }
                inline double& operator[](int i) { return data()[i]; }
            };

            inline void two_sum(double a, double b, double& x, double& y)
            {
                x = a + b;
                double bv = x - a;
                double av = x - bv;
                y = (a - av) + (b - bv);
            }

            inline void two_prod(double a, double b, double& x, double& y)
            {
                x = a * b;
                y = std::fma(a, b, -x);
            }

            inline void diff(double a, double b, Exp& r)
            {
                r.clear();
                double x, y;
                two_sum(a, -b, x, y);
                if (y != 0.0) r.push_back(y);
                if (x != 0.0) r.push_back(x);
            }

            inline void grow(const Exp& e, double b, Exp& r)
            {
                r.clear();
                double q = b;
                for (int i = 0; i < e.size; ++i)
                {
                    double ei = e[i];
                    double qn, h;
                    two_sum(q, ei, qn, h);
                    q = qn;
                    if (h != 0.0) r.push_back(h);
                }
                if (q != 0.0) r.push_back(q);
            }

            inline void add(const Exp& a, const Exp& b, Exp& r)
            {
                r = a;
                Exp temp;
                for (int i = 0; i < b.size; ++i)
                {
                    grow(r, b[i], temp);
                    r = temp;
                }
            }

            inline void neg(const Exp& a, Exp& r)
            {
                r.clear();
                for (int i = 0; i < a.size; ++i) r.push_back(-a[i]);
            }

            inline void sub(const Exp& a, const Exp& b, Exp& r)
            {
                Exp nega_b;
                neg(b, nega_b);
                add(a, nega_b, r);
            }

            inline void scale(const Exp& e, double b, Exp& r)
            {
                r.clear();
                if (e.empty() || b == 0.0) return;
                double q, hh;
                two_prod(e[0], b, q, hh);
                if (hh != 0.0) r.push_back(hh);
                for (int i = 1; i < e.size; ++i)
                {
                    double p1, p0, sum;
                    two_prod(e[i], b, p1, p0);
                    two_sum(q, p0, sum, hh);
                    if (hh != 0.0) r.push_back(hh);
                    two_sum(p1, sum, q, hh);
                    if (hh != 0.0) r.push_back(hh);
                }
                if (q != 0.0) r.push_back(q);
            }

            inline void mul(const Exp& a, const Exp& b, Exp& r)
            {
                r.clear();
                if (a.empty() || b.empty()) return;
                const Exp& big = a.size >= b.size ? a : b;
                const Exp& small = a.size >= b.size ? b : a;
                Exp scaled_big, added_res;
                for (int i = 0; i < small.size; ++i)
                {
                    scale(big, small[i], scaled_big);
                    add(r, scaled_big, added_res);
                    r = added_res;
                }
            }

            inline int sign(const Exp& e)
            {
                if (e.empty()) return 0;
                return e.back() > 0 ? 1 : -1;
            }
        }

        constexpr double EPS = 1.1102230246251565e-16;
        constexpr double ERR_ORIENT = (3.0 + 16.0 * EPS) * EPS;
        constexpr double ERR_INCIRCLE = (10.0 + 96.0 * EPS) * EPS;

#if defined(_MSC_VER)
#define ROBUST_NOINLINE __declspec(noinline)
#else
#define ROBUST_NOINLINE __attribute__((noinline))
#endif

        ROBUST_NOINLINE double orient2d_exact(const Point& a, const Point& b, const Point& c)
        {
            using namespace exact;
            Exp acx, acy, bcx, bcy, m1, m2, s;
            diff(a.x, c.x, acx); diff(a.y, c.y, acy);
            diff(b.x, c.x, bcx); diff(b.y, c.y, bcy);
            mul(acx, bcy, m1);
            mul(acy, bcx, m2);
            sub(m1, m2, s);
            return static_cast<double>(sign(s));
        }

        inline double orient2d(const Point& a, const Point& b, const Point& c)
        {
            const double detleft = (a.x - c.x) * (b.y - c.y);
            const double detright = (a.y - c.y) * (b.x - c.x);
            const double det = detleft - detright;
            const double bound = ERR_ORIENT * (std::fabs(detleft) + std::fabs(detright));
            if (det > bound || -det > bound) return det;
            return orient2d_exact(a, b, c);
        }

        ROBUST_NOINLINE double incircle_exact(const Point& a, const Point& b, const Point& c, const Point& d)
        {
            using namespace exact;
            Exp eadx, eady, ebdx, ebdy, ecdx, ecdy;
            diff(a.x, d.x, eadx); diff(a.y, d.y, eady);
            diff(b.x, d.x, ebdx); diff(b.y, d.y, ebdy);
            diff(c.x, d.x, ecdx); diff(c.y, d.y, ecdy);
            Exp m1, m2, ea, eb, ec;
            mul(eadx, eadx, m1); mul(eady, eady, m2); add(m1, m2, ea);
            mul(ebdx, ebdx, m1); mul(ebdy, ebdy, m2); add(m1, m2, eb);
            mul(ecdx, ecdx, m1); mul(ecdy, ecdy, m2); add(m1, m2, ec);
            Exp sub1, sub2, sub3, mul1, mul2, mul3;
            mul(ebdx, ecdy, m1); mul(ecdx, ebdy, m2); sub(m1, m2, sub1);
            mul(ea, sub1, mul1);
            mul(ecdx, eady, m1); mul(eadx, ecdy, m2); sub(m1, m2, sub2);
            mul(eb, sub2, mul2);
            mul(eadx, ebdy, m1); mul(ebdx, eady, m2); sub(m1, m2, sub3);
            mul(ec, sub3, mul3);
            Exp add1, total_sum;
            add(mul1, mul2, add1);
            add(add1, mul3, total_sum);
            return static_cast<double>(sign(total_sum));
        }

        inline double incircle(const Point& a, const Point& b, const Point& c, const Point& d)
        {
            const double adx = a.x - d.x, ady = a.y - d.y;
            const double bdx = b.x - d.x, bdy = b.y - d.y;
            const double cdx = c.x - d.x, cdy = c.y - d.y;
            const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
            const double cdxady = cdx * ady, adxcdy = adx * cdy;
            const double adxbdy = adx * bdy, bdxady = bdx * ady;
            const double alift = adx * adx + ady * ady;
            const double blift = bdx * bdx + bdy * bdy;
            const double clift = cdx * cdx + cdy * cdy;
            const double det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) +
                clift * (adxbdy - bdxady);
            const double permanent = (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * alift +
                (std::fabs(cdxady) + std::fabs(adxcdy)) * blift +
                (std::fabs(adxbdy) + std::fabs(bdxady)) * clift;
            const double bound = ERR_INCIRCLE * permanent;
            if (det > bound || -det > bound) return det;
            return incircle_exact(a, b, c, d);
        }
    }

    namespace detail
    {
        struct Vtx : Point { size_t id; };

        inline bool lessF(const Point& a, const Point& b, int frame)
        {
            if (frame == 0)
                return a.x < b.x || (a.x == b.x && a.y < b.y);
            return a.y < b.y || (a.y == b.y && a.x > b.x);
        }

        inline int rot(int e) { return (e & ~3) | ((e + 1) & 3); }
        inline int invrot(int e) { return (e & ~3) | ((e + 3) & 3); }
        inline int sym(int e) { return e ^ 2; }

        struct QE
        {
            int next[4];
            int vtx[2];
        };

        class Builder
        {
        public:
            static constexpr int CHUNK_SIZE = 65536;
            static constexpr int MAX_POOLS = 512;
            static constexpr int MAX_PARALLEL_DEPTH = 3;
            static constexpr int PARALLEL_MIN_SIZE = 524288;
            static constexpr int MERGE_PARALLEL_DEPTH = 1;

            Builder(std::vector<Vtx>& v, bool dwyer, ThreadPool* pool)
                : V(v), dwyer_(dwyer), pool_(pool)
            {
                qe.resize(5 * static_cast<int>(V.size()) + 64 * 1024);
                for (auto& p : pools_)
                {
                    p.next_slot = 0;
                    p.block_end = 0;
                    p.free_head = -1;
                }
            }

            void run()
            {
                build(0, static_cast<int>(V.size()), 0, 0, true, 0);
            }

            void extract(std::vector<std::array<size_t, 3>>& out) const
            {
                out.reserve(2 * V.size());
                const int sz = g_next.load(std::memory_order_acquire);
                for (int k = 0; k < sz; ++k)
                {
                    if (qe[k].vtx[0] < 0) continue;
                    for (int r = 0; r < 4; r += 2)
                    {
                        const int e = 4 * k + r;
                        const int f = lnext(e), g = lnext(f);
                        if (lnext(g) != e) continue;
                        if (e > f || e > g) continue;
                        const int a = org(e), b = org(f), c = org(g);
                        if (robust::orient2d(V[a], V[b], V[c]) > 0)
                            out.push_back({ V[a].id, V[b].id, V[c].id });
                    }
                }
            }

        private:
            struct alignas(64) Pool
            {
                int next_slot = 0;
                int block_end = 0;
                int free_head = -1;
                char pad[64 - 3 * sizeof(int)];
            };
            static_assert(sizeof(Pool) == 64, "Pool must fit exactly one cache line");

            std::vector<Vtx>& V;
            std::vector<QE> qe;
            std::array<Pool, MAX_POOLS> pools_{};
            std::atomic<int> g_next{ 0 };
            bool dwyer_;
            ThreadPool* pool_ = nullptr;

            Pool& my_pool()
            {
                const int tid = get_thread_id();
                if (tid < 0 || tid >= MAX_POOLS)
                {
                    std::fprintf(stderr, "Builder: thread id %d >= MAX_POOLS\n", tid);
                    std::abort();
                }
                return pools_[tid];
            }

            int onext(int e) const { return qe[e >> 2].next[e & 3]; }
            int org(int e) const { return qe[e >> 2].vtx[(e >> 1) & 1]; }
            int dest(int e) const { return org(sym(e)); }
            void setOrg(int e, int v) { qe[e >> 2].vtx[(e >> 1) & 1] = v; }

            int oprev(int e) const { return rot(onext(rot(e))); }
            int lnext(int e) const { return rot(onext(invrot(e))); }
            int rprev(int e) const { return onext(sym(e)); }

            int makeEdge()
            {
                Pool& p = my_pool();
                int k;
                if (p.free_head >= 0)
                {
                    k = p.free_head;
                    p.free_head = qe[k].next[0];
                }
                else
                {
                    if (p.next_slot >= p.block_end)
                    {
                        const int s = g_next.fetch_add(CHUNK_SIZE, std::memory_order_acq_rel);
                        if (s + CHUNK_SIZE > static_cast<int>(qe.size()))
                        {
                            std::fprintf(stderr, "Builder: qe overflow (s=%d cap=%zu)\n",
                                s, qe.size());
                            std::abort();
                        }
                        for (int i = 0; i < CHUNK_SIZE; ++i)
                        {
                            qe[s + i].vtx[0] = -1;
                            qe[s + i].vtx[1] = -1;
                        }
                        p.next_slot = s;
                        p.block_end = s + CHUNK_SIZE;
                    }
                    k = p.next_slot++;
                }
                QE& q = qe[k];
                const int e = 4 * k;
                q.next[0] = e;
                q.next[1] = e + 3;
                q.next[2] = e + 2;
                q.next[3] = e + 1;
                q.vtx[0] = q.vtx[1] = -1;
                return e;
            }

            void deleteEdge(int e)
            {
                splice(e, oprev(e));
                splice(sym(e), oprev(sym(e)));
                const int k = e >> 2;
                Pool& p = my_pool();
                qe[k].vtx[0] = -1;
                qe[k].next[0] = p.free_head;
                p.free_head = k;
            }

            void splice(int a, int b)
            {
                const int alpha = rot(onext(a)), beta = rot(onext(b));
                std::swap(qe[a >> 2].next[a & 3], qe[b >> 2].next[b & 3]);
                std::swap(qe[alpha >> 2].next[alpha & 3],
                    qe[beta >> 2].next[beta & 3]);
            }

            int connect(int a, int b)
            {
                const int e = makeEdge();
                setOrg(e, dest(a));
                setOrg(sym(e), org(b));
                splice(e, lnext(a));
                splice(sym(e), b);
                return e;
            }

            bool ccw(int a, int b, int c) const
            {
                return robust::orient2d(V[a], V[b], V[c]) > 0;
            }
            bool rightOf(int x, int e) const { return ccw(x, dest(e), org(e)); }
            bool leftOf(int x, int e) const { return ccw(x, org(e), dest(e)); }

            bool inCircle(int a, int b, int c, int d) const
            {
                if (a == d || b == d || c == d || a == b || a == c || b == c)
                    return false;
                return robust::incircle(V[a], V[b], V[c], V[d]) > 0;
            }

            std::pair<int, int> leaf(int lo, int hi, int frame)
            {
                std::sort(V.begin() + lo, V.begin() + hi,
                    [frame](const Vtx& a, const Vtx& b) { return lessF(a, b, frame); });
                const int n = hi - lo;
                const int s1 = lo, s2 = lo + 1;
                if (n == 2)
                {
                    const int a = makeEdge();
                    setOrg(a, s1);
                    setOrg(sym(a), s2);
                    return { a, sym(a) };
                }
                const int s3 = lo + 2;
                const int a = makeEdge(), b = makeEdge();
                splice(sym(a), b);
                setOrg(a, s1);
                setOrg(sym(a), s2);
                setOrg(b, s2);
                setOrg(sym(b), s3);
                const double o = robust::orient2d(V[s1], V[s2], V[s3]);
                if (o > 0) { connect(b, a); return { a, sym(b) }; }
                if (o < 0) { const int c = connect(b, a); return { sym(c), c }; }
                return { a, sym(b) };
            }

            std::pair<int, bool> merge_scan_left(int basel)
            {
                const int bo = org(basel);
                const int bd = dest(basel);
                auto validE = [&](int e) { return ccw(dest(e), bd, bo); };

                int lc = onext(sym(basel));
                bool lv = validE(lc);
                if (lv)
                {
                    bool deleted = false;
                    while (inCircle(bd, bo, dest(lc), dest(onext(lc))))
                    {
                        const int t = onext(lc);
                        deleteEdge(lc);
                        lc = t;
                        deleted = true;
                    }
                    if (deleted) lv = validE(lc);
                }
                return { lc, lv };
            }

            std::pair<int, bool> merge_scan_right(int basel)
            {
                const int bo = org(basel);
                const int bd = dest(basel);
                auto validE = [&](int e) { return ccw(dest(e), bd, bo); };

                int rc = oprev(basel);
                bool rv = validE(rc);
                if (rv)
                {
                    bool deleted = false;
                    while (inCircle(bd, bo, dest(rc), dest(oprev(rc))))
                    {
                        const int t = oprev(rc);
                        deleteEdge(rc);
                        rc = t;
                        deleted = true;
                    }
                    if (deleted) rv = validE(rc);
                }
                return { rc, rv };
            }

            std::pair<int, int> merge(int ldo, int ldi, int rdi, int rdo, int depth)
            {
                for (;;)
                {
                    if (leftOf(org(rdi), ldi)) ldi = lnext(ldi);
                    else if (rightOf(org(ldi), rdi)) rdi = rprev(rdi);
                    else break;
                }
                int basel = connect(sym(rdi), ldi);
                if (org(ldi) == org(ldo)) ldo = sym(basel);
                if (org(rdi) == org(rdo)) rdo = basel;

                const bool parallel_merge =
                    pool_ != nullptr && depth < MERGE_PARALLEL_DEPTH;

                for (;;)
                {
                    int lc, rc;
                    bool lv, rv;

                    if (parallel_merge)
                    {
                        auto fut = pool_->enqueue(
                            [this, basel]() { return merge_scan_left(basel); });
                        auto r = merge_scan_right(basel);
                        rc = r.first;
                        rv = r.second;
                        auto l = pool_->wait_for(fut);
                        lc = l.first;
                        lv = l.second;
                    }
                    else
                    {
                        auto l = merge_scan_left(basel);
                        lc = l.first;
                        lv = l.second;
                        auto r = merge_scan_right(basel);
                        rc = r.first;
                        rv = r.second;
                    }

                    if (!lv && !rv) break;

                    const int bo = org(basel);
                    const int bd = dest(basel);
                    if (!lv || (rv && inCircle(dest(lc), bd, bo, dest(rc))))
                        basel = connect(rc, sym(basel));
                    else
                        basel = connect(sym(basel), sym(lc));
                }
                return { ldo, rdo };
            }

            int hullWalk(int e, int frame, bool wantMax) const
            {
                auto better = [&](int u, int v)
                    {
                        return wantMax ? lessF(V[v], V[u], frame)
                            : lessF(V[u], V[v], frame);
                    };
                for (;;)
                {
                    const int nx = rprev(e);
                    const int pv = sym(oprev(e));
                    if (better(org(nx), org(e))) e = nx;
                    else if (better(org(pv), org(e))) e = pv;
                    else return e;
                }
            }

            std::pair<int, int> build(int lo, int hi, int axis, int outAxis,
                bool presorted, int depth)
            {
                const int n = hi - lo;
                if (n <= 3) return leaf(lo, hi, outAxis);

                const int mid = lo + n / 2;
                if (dwyer_ && !presorted)
                {
                    if (axis == 0)
                    {
                        std::nth_element(V.begin() + lo, V.begin() + mid, V.begin() + hi,
                            [](const Vtx& a, const Vtx& b) {
                                return a.x < b.x || (a.x == b.x && a.y < b.y);
                            });
                    }
                    else
                    {
                        std::nth_element(V.begin() + lo, V.begin() + mid, V.begin() + hi,
                            [](const Vtx& a, const Vtx& b) {
                                return a.y < b.y || (a.y == b.y && a.x > b.x);
                            });
                    }
                }
                const int childAxis = dwyer_ ? 1 - axis : axis;

                std::pair<int, int> L, R;
                const bool can_parallel =
                    pool_ != nullptr &&
                    depth < MAX_PARALLEL_DEPTH &&
                    n >= PARALLEL_MIN_SIZE;

                if (can_parallel)
                {
                    auto fut = pool_->enqueue(
                        [this, lo, mid, childAxis, axis, depth]() {
                            return build(lo, mid, childAxis, axis, false, depth + 1);
                        });
                    R = build(mid, hi, childAxis, axis, false, depth + 1);
                    L = pool_->wait_for(fut);
                }
                else
                {
                    L = build(lo, mid, childAxis, axis, false, depth + 1);
                    R = build(mid, hi, childAxis, axis, false, depth + 1);
                }

                auto res = merge(L.first, L.second, R.first, R.second, depth);

                if (outAxis != axis)
                {
                    const int le = hullWalk(res.first, outAxis, false);
                    const int reCcw = hullWalk(onext(res.second), outAxis, true);
                    res.first = le;
                    res.second = oprev(reCcw);
                }
                return res;
            }
        };
    }

    wrapper::Triangulation delaunay_triangulate(Point* pts, size_t sz, bool dwyer)
    {
        const int n = static_cast<int>(sz);

        wrapper::Triangulation res;
        res.rep.assign(n, 0);

        std::vector<detail::Vtx> V;
        V.reserve(n);
        for (int i = 0; i < n; ++i)
        {
            detail::Vtx v;
            v.x = pts[i].x;
            v.y = pts[i].y;
            v.id = static_cast<size_t>(i);
            V.push_back(v);
        }
        std::iota(res.rep.begin(), res.rep.end(), static_cast<size_t>(0));
        std::sort(V.begin(), V.end(), [](const detail::Vtx& a, const detail::Vtx& b)
            {
                if (a.x != b.x) return a.x < b.x;
                if (a.y != b.y) return a.y < b.y;
                return a.id < b.id;
            });
        {
            size_t m = 0;
            for (size_t i = 0; i < V.size(); ++i)
            {
                if (m > 0 && V[m - 1].x == V[i].x && V[m - 1].y == V[i].y)
                {
                    res.rep[V[i].id] = V[m - 1].id;
                    continue;
                }
                if (m != i) V[m] = V[i];
                ++m;
            }
            V.erase(V.begin() + m, V.end());
        }
        res.unique_points = static_cast<int>(V.size());

        if (V.size() >= 3)
        {
            unsigned nthreads = std::thread::hardware_concurrency();
            if (nthreads == 0) nthreads = 1;
            if (nthreads > 64) nthreads = 64;

            if (nthreads > 1 && V.size() >= 200000)
            {
                ThreadPool pool(nthreads);
                detail::Builder b(V, dwyer, &pool);
                b.run();
                b.extract(res.triangles);
            }
            else
            {
                detail::Builder b(V, dwyer, nullptr);
                b.run();
                b.extract(res.triangles);
            }
        }
        return res;
    }
}
