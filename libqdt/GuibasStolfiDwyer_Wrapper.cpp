#include <algorithm>
#include <cstdio>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "qdt.h"
#include "wrapper.hxx"

namespace wrapper
{
    using Point = qdt_vec2d;
    Triangulation delaunay_triangulate(Point* pts, size_t sz, bool dwyer);
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

                ~Exp()
                {
                    if (heap_ptr)
                        delete[] heap_ptr;
                }

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
                    {
                        std::copy(other.local_buf, other.local_buf + size, local_buf);
                    }
                    other.heap_ptr = nullptr;
                    other.size = 0;
                    other.capacity = LOCAL_CAP;
                }

                Exp& operator=(const Exp& other)
                {
                    if (this == &other)
                        return *this;
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
                        if (heap_ptr)
                        {
                            delete[] heap_ptr;
                            heap_ptr = nullptr;
                        }
                        capacity = LOCAL_CAP;
                        std::copy(other.local_buf, other.local_buf + other.size, local_buf);
                    }
                    size = other.size;
                    return *this;
                }

                Exp& operator=(Exp&& other) noexcept
                {
                    if (this == &other)
                        return *this;
                    if (heap_ptr)
                        delete[] heap_ptr;
                    size = other.size;
                    capacity = other.capacity;
                    heap_ptr = other.heap_ptr;
                    if (!heap_ptr)
                    {
                        std::copy(other.local_buf, other.local_buf + size, local_buf);
                    }
                    other.heap_ptr = nullptr;
                    other.size = 0;
                    other.capacity = LOCAL_CAP;
                    return *this;
                }

                inline void clear()
                {
                    size = 0;
                }

                inline double* data()
                {
                    return heap_ptr ? heap_ptr : local_buf;
                }

                inline const double* data() const
                {
                    return heap_ptr ? heap_ptr : local_buf;
                }

                inline void push_back(double val)
                {
                    if (size >= capacity)
                    {
                        int new_cap = capacity * 2;
                        double* new_heap = new double[new_cap];
                        std::copy(data(), data() + size, new_heap);
                        if (heap_ptr)
                            delete[] heap_ptr;
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
                if (y != 0.0)
                    r.push_back(y);
                if (x != 0.0)
                    r.push_back(x);
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
                    if (h != 0.0)
                        r.push_back(h);
                }
                if (q != 0.0)
                    r.push_back(q);
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
                for (int i = 0; i < a.size; ++i)
                    r.push_back(-a[i]);
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
                if (e.empty() || b == 0.0)
                    return;
                double q, hh;
                two_prod(e[0], b, q, hh);
                if (hh != 0.0)
                    r.push_back(hh);
                for (int i = 1; i < e.size; ++i)
                {
                    double p1, p0, sum;
                    two_prod(e[i], b, p1, p0);
                    two_sum(q, p0, sum, hh);
                    if (hh != 0.0)
                        r.push_back(hh);
                    two_sum(p1, sum, q, hh);
                    if (hh != 0.0)
                        r.push_back(hh);
                }
                if (q != 0.0)
                    r.push_back(q);
            }

            inline void mul(const Exp& a, const Exp& b, Exp& r)
            {
                r.clear();
                if (a.empty() || b.empty())
                    return;
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
                if (e.empty())
                    return 0;
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
            if (det > bound || -det > bound)
                return det;
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
            const double det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
            const double permanent = (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * alift +
                (std::fabs(cdxady) + std::fabs(adxcdy)) * blift +
                (std::fabs(adxbdy) + std::fabs(bdxady)) * clift;
            const double bound = ERR_INCIRCLE * permanent;
            if (det > bound || -det > bound)
                return det;
            return incircle_exact(a, b, c, d);
        }
    }
    namespace detail
    {
        struct Vtx : Point
        {
            size_t id;
        };

        inline bool lessF(const Point& a, const Point& b, int frame)
        {
            if (frame == 0)
                return a.x < b.x || (a.x == b.x && a.y < b.y);
            return a.y < b.y || (a.y == b.y && a.x > b.x);
        }

        inline int rot(int e)
        {
            return (e & ~3) | ((e + 1) & 3);
        }
        inline int invrot(int e)
        {
            return (e & ~3) | ((e + 3) & 3);
        }
        inline int sym(int e)
        {
            return e ^ 2;
        }

        struct QE
        {
            int next[4];
            int vtx[2];
        };

        class Builder
        {
        public:
            Builder(std::vector<Vtx>& v, bool dwyer) : V(v), dwyer_(dwyer)
            {
                qe.reserve(3 * V.size() + 16);
            }

            void run()
            {
                build(0, static_cast<int>(V.size()), 0, 0, true);
            }

            void extract(std::vector<std::array<size_t, 3>>& out) const
            {
                out.reserve(2 * V.size());
                for (size_t k = 0; k < qe.size(); ++k)
                {
                    if (qe[k].vtx[0] < 0)
                        continue;
                    for (int r = 0; r < 4; r += 2)
                    {
                        const int e = static_cast<int>(4 * k) + r;
                        const int f = lnext(e), g = lnext(f);
                        if (lnext(g) != e)
                            continue;
                        if (e > f || e > g)
                            continue;
                        const int a = org(e), b = org(f), c = org(g);
                        if (robust::orient2d(V[a], V[b], V[c]) > 0)
                            out.push_back({ V[a].id, V[b].id, V[c].id });
                    }
                }
            }

        private:
            std::vector<Vtx>& V;
            std::vector<QE> qe;
            int freeHead = -1;
            bool dwyer_;

            int onext(int e) const
            {
                return qe[e >> 2].next[e & 3];
            }
            int org(int e) const
            {
                return qe[e >> 2].vtx[(e >> 1) & 1];
            }
            int dest(int e) const
            {
                return org(sym(e));
            }
            void setOrg(int e, int v)
            {
                qe[e >> 2].vtx[(e >> 1) & 1] = v;
            }

            int oprev(int e) const
            {
                return rot(onext(rot(e)));
            }
            int lnext(int e) const
            {
                return rot(onext(invrot(e)));
            }
            int rprev(int e) const
            {
                return onext(sym(e));
            }

            int makeEdge()
            {
                int k;
                if (freeHead >= 0)
                {
                    k = freeHead;
                    freeHead = qe[k].next[0];
                }
                else
                {
                    k = static_cast<int>(qe.size());
                    qe.emplace_back();
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

            void splice(int a, int b)
            {
                const int alpha = rot(onext(a)), beta = rot(onext(b));
                std::swap(qe[a >> 2].next[a & 3], qe[b >> 2].next[b & 3]);
                std::swap(qe[alpha >> 2].next[alpha & 3], qe[beta >> 2].next[beta & 3]);
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

            void deleteEdge(int e) {
                splice(e, oprev(e));
                splice(sym(e), oprev(sym(e)));
                const int k = e >> 2;
                qe[k].vtx[0] = -1;
                qe[k].next[0] = freeHead;
                freeHead = k;
            }

            bool ccw(int a, int b, int c) const
            {
                return robust::orient2d(V[a], V[b], V[c]) > 0;
            }
            bool rightOf(int x, int e) const
            {
                return ccw(x, dest(e), org(e));
            }
            bool leftOf(int x, int e) const
            {
                return ccw(x, org(e), dest(e));
            }
            bool inCircle(int a, int b, int c, int d) const
            {
                if (a == d || b == d || c == d || a == b || a == c || b == c)
                    return false;
                return robust::incircle(V[a], V[b], V[c], V[d]) > 0;
            }
            bool valid(int e, int basel) const
            {
                return rightOf(dest(e), basel);
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
                if (o > 0)
                {
                    connect(b, a);
                    return { a, sym(b) };
                }
                if (o < 0)
                {
                    const int c = connect(b, a);
                    return { sym(c), c };
                }
                return { a, sym(b) };
            }

            std::pair<int, int> merge(int ldo, int ldi, int rdi, int rdo)
            {
                for (;;)
                {
                    if (leftOf(org(rdi), ldi))
                        ldi = lnext(ldi);
                    else if (rightOf(org(ldi), rdi))
                        rdi = rprev(rdi);
                    else
                        break;
                }
                int basel = connect(sym(rdi), ldi);
                if (org(ldi) == org(ldo))
                    ldo = sym(basel);
                if (org(rdi) == org(rdo))
                    rdo = basel;

                for (;;)
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
                    if (!lv && !rv)
                        break;
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
                        return wantMax ? lessF(V[v], V[u], frame) : lessF(V[u], V[v], frame);
                    };
                for (;;)
                {
                    const int nx = rprev(e);
                    const int pv = sym(oprev(e));
                    if (better(org(nx), org(e)))
                        e = nx;
                    else if (better(org(pv), org(e)))
                        e = pv;
                    else
                        return e;
                }
            }

            std::pair<int, int> build(int lo, int hi, int axis, int outAxis, bool presorted)
            {
                const int n = hi - lo;
                if (n <= 3)
                    return leaf(lo, hi, outAxis);

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
                const auto L = build(lo, mid, childAxis, axis, false);
                const auto R = build(mid, hi, childAxis, axis, false);
                auto res = merge(L.first, L.second, R.first, R.second);

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
    void Triangulation::ensure_mesh_()
    {
        if (mesh_built_) return;

        std::unordered_map<size_t, int> orig_to_local;
        mesh_vertices_.clear();
        mesh_halfedges_.clear();
        free_halfedges_.clear();
        last_face_he_ = 0;

        for (const auto& t : triangles)
        {
            for (int i = 0; i < 3; ++i)
            {
                const size_t orig = t[i];
                if (orig_to_local.find(orig) == orig_to_local.end())
                {
                    const int local = static_cast<int>(mesh_vertices_.size());
                    orig_to_local[orig] = local;

                    MeshVertex v;
                    v.p = input_points[orig];
                    v.orig_id = orig;
                    v.halfedge = -1;
                    mesh_vertices_.push_back(v);
                }
            }
        }

        std::unordered_map<uint64_t, int> edge_map;

        for (const auto& t : triangles)
        {
            int v[3];
            for (int i = 0; i < 3; ++i)
                v[i] = orig_to_local[t[i]];

            int h[3];
            for (int i = 0; i < 3; ++i)
            {
                h[i] = static_cast<int>(mesh_halfedges_.size());
                mesh_halfedges_.push_back({ v[i], -1, -1 });
            }

            for (int i = 0; i < 3; ++i)
                mesh_halfedges_[h[i]].next = h[(i + 1) % 3];

            for (int i = 0; i < 3; ++i)
            {
                const int a = v[i];
                const int b = v[(i + 1) % 3];

                const uint64_t key_fwd =
                    (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
                    static_cast<uint32_t>(b);
                const uint64_t key_rev =
                    (static_cast<uint64_t>(static_cast<uint32_t>(b)) << 32) |
                    static_cast<uint32_t>(a);

                auto it = edge_map.find(key_rev);
                if (it != edge_map.end())
                {
                    const int other = it->second;
                    mesh_halfedges_[h[i]].twin = other;
                    mesh_halfedges_[other].twin = h[i];
                    edge_map.erase(it);
                }
                else
                {
                    edge_map[key_fwd] = h[i];
                }
            }

            for (int i = 0; i < 3; ++i)
            {
                if (mesh_vertices_[v[i]].halfedge < 0)
                    mesh_vertices_[v[i]].halfedge = h[i];
            }
        }

        mesh_built_ = true;
    }

    int Triangulation::alloc_halfedge_()
    {
        if (!free_halfedges_.empty())
        {
            const int h = free_halfedges_.back();
            free_halfedges_.pop_back();
            return h;
        }
        const int h = static_cast<int>(mesh_halfedges_.size());
        mesh_halfedges_.push_back({});
        return h;
    }

    void Triangulation::free_halfedge_(int h)
    {
        mesh_halfedges_[h].origin = -1;
        mesh_halfedges_[h].twin = -1;
        mesh_halfedges_[h].next = -1;
        free_halfedges_.push_back(h);
    }

    // 定位包含 p 的三角形，返回该三角形的一条半边。
    // 若 p 位于凸包之外，返回 -1，并通过 hull_he 给出 walk 越过的那条凸包半边
    // （p 严格位于该边外侧，即该边对 p 可见）。
    // 返回 -1 且 hull_he == -1 表示定位失败（理论上不应发生）。
    int Triangulation::locate_(const qdt_vec2d& p, int& hull_he)
    {
        hull_he = -1;
        if (mesh_halfedges_.empty())
            return -1;

        int h = last_face_he_;
        if (h < 0 || h >= static_cast<int>(mesh_halfedges_.size()) ||
            mesh_halfedges_[h].origin < 0)
        {
            h = 0;
            while (h < static_cast<int>(mesh_halfedges_.size()) &&
                mesh_halfedges_[h].origin < 0)
                ++h;
            if (h >= static_cast<int>(mesh_halfedges_.size()))
                return -1;
        }

        for (int step = 0; step < 1000000; ++step)
        {
            const int h0 = h;
            const int h1 = mesh_halfedges_[h0].next;
            const int h2 = mesh_halfedges_[h1].next;
            if (mesh_halfedges_[h2].next != h0)
                break;

            const int v0 = mesh_halfedges_[h0].origin;
            const int v1 = mesh_halfedges_[h1].origin;
            const int v2 = mesh_halfedges_[h2].origin;

            const qdt_vec2d& A = mesh_vertices_[v0].p;
            const qdt_vec2d& B = mesh_vertices_[v1].p;
            const qdt_vec2d& C = mesh_vertices_[v2].p;

            const double o0 = robust::orient2d(A, B, p);
            const double o1 = robust::orient2d(B, C, p);
            const double o2 = robust::orient2d(C, A, p);

            if (o0 >= 0 && o1 >= 0 && o2 >= 0)
            {
                last_face_he_ = h0;
                return h0;
            }

            int cross;
            if (o0 < 0)      cross = h0;
            else if (o1 < 0) cross = h1;
            else             cross = h2;

            const int tw = mesh_halfedges_[cross].twin;
            if (tw < 0)
            {
                hull_he = cross;
                return -1;
            }

            h = tw;
        }

        return -1;
    }

    void Triangulation::rebuild_output_()
    {
        triangles.clear();
        std::unordered_set<int> visited;

        for (int h = 0; h < static_cast<int>(mesh_halfedges_.size()); ++h)
        {
            if (mesh_halfedges_[h].origin < 0)
                continue;

            const int h1 = mesh_halfedges_[h].next;
            const int h2 = mesh_halfedges_[h1].next;
            if (h1 < 0 || h2 < 0)
                continue;
            if (mesh_halfedges_[h2].next != h)
                continue;

            const int canon = std::min({ h, h1, h2 });
            if (visited.count(canon))
                continue;
            visited.insert(canon);

            const int v0 = mesh_halfedges_[h].origin;
            const int v1 = mesh_halfedges_[h1].origin;
            const int v2 = mesh_halfedges_[h2].origin;

            triangles.push_back({
                mesh_vertices_[v0].orig_id,
                mesh_vertices_[v1].orig_id,
                mesh_vertices_[v2].orig_id
                });
        }
    }

    // 兜底路径：把点追加到 input_points 后整体重建。
    // 用于网格为空（有效点不足 3 个或全部共线）以及定位失败等罕见情况，保证结果始终正确。
    void Triangulation::insert_by_rebuild_(const qdt_vec2d& p)
    {
        input_points.push_back(p);
        Triangulation fresh = delaunay_triangulate(
            input_points.data(), input_points.size(), true);
        *this = std::move(fresh);
    }

    void Triangulation::insert(const qdt_vec2d& p)
    {
        ensure_mesh_();

        for (const auto& v : mesh_vertices_)
        {
            if (v.p.x == p.x && v.p.y == p.y)
            {
                input_points.push_back(p);
                rep.push_back(v.orig_id);
                return;
            }
        }

        // 还没有三角形（点数不足 3 个或全部共线）：直接整体重建
        if (mesh_halfedges_.empty())
        {
            insert_by_rebuild_(p);
            return;
        }

        int hull_he = -1;
        const int h_start = locate_(p, hull_he);
        if (h_start < 0 && hull_he < 0)
        {
            insert_by_rebuild_(p);
            return;
        }
        const bool outside = (h_start < 0);   // p 在当前凸包之外

        const size_t new_id = input_points.size();
        input_points.push_back(p);
        rep.push_back(new_id);
        ++unique_points;

        const int new_v = static_cast<int>(mesh_vertices_.size());
        {
            MeshVertex v;
            v.p = p;
            v.orig_id = new_id;
            v.halfedge = -1;
            mesh_vertices_.push_back(v);
        }

        auto he_dest = [&](int e) {
            return mesh_halfedges_[mesh_halfedges_[e].next].origin;
            };

        // ---- p 在凸包外：收集所有对 p 可见的凸包边（它们在凸包上是连续的一段） ----
        std::vector<int> visible_hull;
        if (outside)
        {
            auto is_visible = [&](int e) {
                const qdt_vec2d& A = mesh_vertices_[mesh_halfedges_[e].origin].p;
                const qdt_vec2d& B = mesh_vertices_[he_dest(e)].p;
                return robust::orient2d(A, B, p) < 0;
                };
            // 凸包上 e 的下一条边（起点 = e 的终点）：绕终点旋转直到碰到无 twin 的半边
            auto next_hull = [&](int e) {
                int x = mesh_halfedges_[e].next;
                while (mesh_halfedges_[x].twin >= 0)
                    x = mesh_halfedges_[mesh_halfedges_[x].twin].next;
                return x;
                };
            // 凸包上 e 的上一条边（终点 = e 的起点）
            auto prev_hull = [&](int e) {
                int x = mesh_halfedges_[mesh_halfedges_[e].next].next;
                while (mesh_halfedges_[x].twin >= 0)
                {
                    const int tw = mesh_halfedges_[x].twin;
                    x = mesh_halfedges_[mesh_halfedges_[tw].next].next;
                }
                return x;
                };

            visible_hull.push_back(hull_he);
            for (int e = next_hull(hull_he); e != hull_he && is_visible(e); e = next_hull(e))
                visible_hull.push_back(e);
            for (int e = prev_hull(hull_he); e != hull_he && is_visible(e); e = prev_hull(e))
                visible_hull.push_back(e);
        }

        // ---- 外接圆包含 p 的三角形构成空洞（BFS） ----
        std::vector<int>        cavity_hes;
        std::unordered_set<int> cavity_set;
        std::queue<int>         q;
        if (!outside)
            q.push(h_start);
        else
            for (int e : visible_hull)
                q.push(e);

        while (!q.empty())
        {
            const int h = q.front(); q.pop();

            const int h1 = mesh_halfedges_[h].next;
            const int h2 = mesh_halfedges_[h1].next;
            const int canon = std::min({ h, h1, h2 });

            if (cavity_set.count(canon))
                continue;

            const int v0 = mesh_halfedges_[h].origin;
            const int v1 = mesh_halfedges_[h1].origin;
            const int v2 = mesh_halfedges_[h2].origin;

            const double ic = robust::incircle(
                mesh_vertices_[v0].p,
                mesh_vertices_[v1].p,
                mesh_vertices_[v2].p,
                p);

            if (ic <= 0)
                continue;

            cavity_set.insert(canon);
            cavity_hes.push_back(h);

            int e = h;
            do
            {
                const int tw = mesh_halfedges_[e].twin;
                if (tw >= 0)
                    q.push(tw);
                e = mesh_halfedges_[e].next;
            } while (e != h);
        }

        // p 在凸包内时，所在三角形必然在空洞中（容错）；凸包外时空洞可以为空
        if (!outside && cavity_hes.empty())
        {
            const int h1 = mesh_halfedges_[h_start].next;
            const int h2 = mesh_halfedges_[h1].next;
            cavity_set.insert(std::min({ h_start, h1, h2 }));
            cavity_hes.push_back(h_start);
        }

        // ---- 收集需要与 p 相连的边，每条边生成一个新三角形 (a, b, p) ----
        struct BoundaryInfo { int a, b, external_twin; };
        std::vector<BoundaryInfo> boundary_info;

        // (1) 空洞边界：空洞三角形的边，若对面没有三角形或对面不在空洞中
        for (int h : cavity_hes)
        {
            int e = h;
            do
            {
                const int tw = mesh_halfedges_[e].twin;
                bool tw_in_cavity = false;
                if (tw >= 0)
                {
                    const int t1 = mesh_halfedges_[tw].next;
                    const int t2 = mesh_halfedges_[t1].next;
                    tw_in_cavity = cavity_set.count(std::min({ tw, t1, t2 })) != 0;
                }

                if (!tw_in_cavity)
                {
                    const int a = mesh_halfedges_[e].origin;
                    const int b = mesh_halfedges_[mesh_halfedges_[e].next].origin;
                    // 凸包边（tw < 0）只有当 p 严格位于其左侧时才能与 p 组成新三角形；
                    // p 在其外侧（边被 p 看到）或恰好落在该边上时都不能生成（后者会退化）
                    if (tw >= 0 ||
                        robust::orient2d(mesh_vertices_[a].p, mesh_vertices_[b].p, p) > 0)
                        boundary_info.push_back({ a, b, tw });
                }

                e = mesh_halfedges_[e].next;
            } while (e != h);
        }

        // (2) 对 p 可见、且所在三角形不在空洞中的凸包边 a->b：
        //     新三角形 (b, a, p)，其边 b->a 与凸包半边 a->b 互为 twin
        for (int e : visible_hull)
        {
            const int e1 = mesh_halfedges_[e].next;
            const int e2 = mesh_halfedges_[e1].next;
            if (cavity_set.count(std::min({ e, e1, e2 })))
                continue;
            const int a = mesh_halfedges_[e].origin;
            const int b = mesh_halfedges_[e1].origin;
            boundary_info.push_back({ b, a, e });
        }

        for (int h : cavity_hes)
        {
            int e = h;
            do
            {
                const int nxt = mesh_halfedges_[e].next;
                free_halfedge_(e);
                e = nxt;
            } while (e != h);
        }

        struct NewTri { int h0, h1, h2; };
        std::vector<NewTri> new_tris;

        for (const auto& bi : boundary_info)
        {
            const int h0 = alloc_halfedge_();
            const int h1 = alloc_halfedge_();
            const int h2 = alloc_halfedge_();

            mesh_halfedges_[h0] = { bi.a, bi.external_twin, h1 };
            mesh_halfedges_[h1] = { bi.b, -1, h2 };
            mesh_halfedges_[h2] = { new_v, -1, h0 };

            if (bi.external_twin >= 0)
                mesh_halfedges_[bi.external_twin].twin = h0;

            new_tris.push_back({ h0, h1, h2 });
        }

        if (new_tris.empty())
        {
            // 不应发生；此时网格已被改动，直接整体重建以保证结果正确
            Triangulation fresh = delaunay_triangulate(
                input_points.data(), input_points.size(), true);
            *this = std::move(fresh);
            return;
        }

        // 新三角形围绕 p 成扇形：h1 (v->p) 与以 v 为起点的另一个新三角形的 h2 (p->v) 互为 twin。
        // 凸包外插入时扇形不闭合，两端留下的 twin = -1 的边就是新的凸包边。
        std::unordered_map<int, int> v_to_newv;
        std::unordered_map<int, int> newv_to_v;

        for (const auto& nt : new_tris)
        {
            const int a = mesh_halfedges_[nt.h0].origin;
            const int b = mesh_halfedges_[nt.h1].origin;
            v_to_newv[b] = nt.h1;
            newv_to_v[a] = nt.h2;
        }

        for (const auto& kv : v_to_newv)
        {
            const int v = kv.first;
            const int h1 = kv.second;
            auto it = newv_to_v.find(v);
            if (it != newv_to_v.end())
            {
                const int h2 = it->second;
                mesh_halfedges_[h1].twin = h2;
                mesh_halfedges_[h2].twin = h1;
            }
        }

        mesh_vertices_[new_v].halfedge = new_tris[0].h2;
        last_face_he_ = new_tris[0].h0;

        rebuild_output_();
    }
    Triangulation delaunay_triangulate(Point* pts, size_t sz, bool dwyer)
    {
        const int n = static_cast<int>(sz);

        Triangulation res;
        res.rep.assign(n, 0);

        std::vector<detail::Vtx> V;
        V.reserve(n);
        for (int i = 0; i < n; ++i) { detail::Vtx v; v.x = pts[i].x; v.y = pts[i].y; v.id = static_cast<size_t>(i); V.push_back(v); }
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
                if (m > 0 && V[m - 1].x == V[i].x && V[m - 1].y == V[i].y) { res.rep[V[i].id] = V[m - 1].id; continue; }
                if (m != i) V[m] = V[i];
                ++m;
            }
            V.erase(V.begin() + m, V.end());
        }
        res.unique_points = static_cast<int>(V.size());

        if (V.size() >= 3)
        {
            detail::Builder b(V, dwyer);
            b.run();
            b.extract(res.triangles);
        }

        res.input_points.resize(sz);
        for (size_t i = 0; i < sz; ++i)
            res.input_points[i] = pts[i];

        return res;
    }
}