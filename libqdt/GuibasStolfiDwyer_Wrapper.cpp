#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <thread>
#include <unordered_map>
#include <vector>

#include "qdt.h"
#include "wrapper.hxx"

namespace wrapper
{
    using Point = qdt_vec2d;
    Triangulation delaunay_triangulate(Point* pts, size_t sz);

    namespace robust
    {
        namespace exact
        {
            struct Exp
            {
                static constexpr int LOCAL_CAP = 8;
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

        // The error bounds of the floating point filters, and the exact
        // arithmetic below, assume that no intermediate product overflows or
        // becomes subnormal. A filter result whose magnitude is below this
        // floor is therefore not trusted, and the exact path is used instead.
        constexpr double FILTER_FLOOR = 1e-250;

        // The exact path first multiplies every coordinate by the same power of
        // two, so that the largest one is about 1. That leaves the sign of the
        // determinant unchanged (orient2d scales by a positive factor, incircle
        // likewise), is exact, and keeps all products far from overflow and
        // underflow, whatever the overall scale of the input. Coordinates whose
        // magnitude is already moderate are left alone.
        inline int normalizing_exponent(const double* v, int n)
        {
            double m = 0.0;
            for (int i = 0; i < n; ++i)
                m = std::max(m, std::fabs(v[i]));
            if (m == 0.0 || !std::isfinite(m)) // all zero, or not a number: leave as is
                return 0;
            const int e = std::ilogb(m);
            return (e >= -64 && e <= 64) ? 0 : e;
        }

        inline Point scaled_point(const Point& p, int e)
        {
            Point r;
            r.x = std::ldexp(p.x, -e);
            r.y = std::ldexp(p.y, -e);
            return r;
        }

#if defined(_MSC_VER)
#define ROBUST_NOINLINE __declspec(noinline)
#else
#define ROBUST_NOINLINE __attribute__((noinline))
#endif
        ROBUST_NOINLINE double orient2d_exact(const Point& a0, const Point& b0, const Point& c0)
        {
            using namespace exact;
            const double raw[6] = { a0.x, a0.y, b0.x, b0.y, c0.x, c0.y };
            const int e = normalizing_exponent(raw, 6);
            const Point a = e ? scaled_point(a0, e) : a0;
            const Point b = e ? scaled_point(b0, e) : b0;
            const Point c = e ? scaled_point(c0, e) : c0;
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
            const double mag = std::fabs(detleft) + std::fabs(detright);
            const double bound = ERR_ORIENT * mag;
            if ((det > bound || -det > bound) && mag > FILTER_FLOOR)
                return det;
            return orient2d_exact(a, b, c);
        }
        ROBUST_NOINLINE double incircle_exact(const Point& a0, const Point& b0, const Point& c0, const Point& d0)
        {
            using namespace exact;
            const double raw[8] = { a0.x, a0.y, b0.x, b0.y, c0.x, c0.y, d0.x, d0.y };
            const int e = normalizing_exponent(raw, 8);
            const Point a = e ? scaled_point(a0, e) : a0;
            const Point b = e ? scaled_point(b0, e) : b0;
            const Point c = e ? scaled_point(c0, e) : c0;
            const Point d = e ? scaled_point(d0, e) : d0;
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
            if ((det > bound || -det > bound) && permanent > FILTER_FLOOR)
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
            Builder(std::vector<Vtx>& v) : V(v) { qe.reserve(3 * V.size() + 16); }
            void run() { build(0, static_cast<int>(V.size()), 0, 0, true); }
            void extract(std::vector<std::array<size_t, 3>>& out) const
            {
                out.reserve(2 * V.size());
                for (size_t k = 0; k < qe.size(); ++k)
                {
                    if (qe[k].vtx[0] < 0) continue;
                    for (int r = 0; r < 4; r += 2)
                    {
                        const int e = static_cast<int>(4 * k) + r;
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
            std::vector<Vtx>& V;
            std::vector<QE> qe;
            int freeHead = -1;

            int onext(int e) const { return qe[e >> 2].next[e & 3]; }
            int org(int e) const { return qe[e >> 2].vtx[(e >> 1) & 1]; }
            int dest(int e) const { return org(sym(e)); }
            void setOrg(int e, int v) { qe[e >> 2].vtx[(e >> 1) & 1] = v; }
            int oprev(int e) const { return rot(onext(rot(e))); }
            int lnext(int e) const { return rot(onext(invrot(e))); }
            int rprev(int e) const { return onext(sym(e)); }

            int makeEdge()
            {
                int k;
                if (freeHead >= 0) { k = freeHead; freeHead = qe[k].next[0]; }
                else { k = static_cast<int>(qe.size()); qe.emplace_back(); }
                QE& q = qe[k];
                const int e = 4 * k;
                q.next[0] = e; q.next[1] = e + 3; q.next[2] = e + 2; q.next[3] = e + 1;
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

            bool ccw(int a, int b, int c) const { return robust::orient2d(V[a], V[b], V[c]) > 0; }
            bool rightOf(int x, int e) const { return ccw(x, dest(e), org(e)); }
            bool leftOf(int x, int e) const { return ccw(x, org(e), dest(e)); }
            bool inCircle(int a, int b, int c, int d) const
            {
                if (a == d || b == d || c == d || a == b || a == c || b == c) return false;
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
                    setOrg(a, s1); setOrg(sym(a), s2);
                    return { a, sym(a) };
                }
                const int s3 = lo + 2;
                const int a = makeEdge(), b = makeEdge();
                splice(sym(a), b);
                setOrg(a, s1); setOrg(sym(a), s2); setOrg(b, s2); setOrg(sym(b), s3);
                const double o = robust::orient2d(V[s1], V[s2], V[s3]);
                if (o > 0) { connect(b, a); return { a, sym(b) }; }
                if (o < 0) { const int c = connect(b, a); return { sym(c), c }; }
                return { a, sym(b) };
            }

            std::pair<int, int> merge(int ldo, int ldi, int rdi, int rdo)
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
                    if (!lv && !rv) break;
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
                    if (better(org(nx), org(e))) e = nx;
                    else if (better(org(pv), org(e))) e = pv;
                    else return e;
                }
            }

            std::pair<int, int> build(int lo, int hi, int axis, int outAxis, bool presorted)
            {
                const int n = hi - lo;
                if (n <= 3) return leaf(lo, hi, outAxis);

                const int mid = lo + n / 2;
                if (!presorted)
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
                const int childAxis = 1 - axis;
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

    size_t Triangulation::add_triangle_to_output_(size_t id0, size_t id1, size_t id2, int h0)
    {
        size_t new_idx = triangles.size();
        triangles.push_back({ id0, id1, id2 });
        triangle_to_he_.push_back(h0);
        return new_idx;
    }

    void Triangulation::remove_triangle_by_index_(size_t idx)
    {
        size_t last_idx = triangles.size() - 1;
        if (idx != last_idx)
        {
            triangles[idx] = triangles[last_idx];
            int moved_he = triangle_to_he_[last_idx];
            triangle_to_he_[idx] = moved_he;

            int cur = moved_he;
            do {
                mesh_halfedges_[cur].triangle_idx = static_cast<int>(idx);
                cur = mesh_halfedges_[cur].next;
            } while (cur != moved_he);
        }
        triangles.pop_back();
        triangle_to_he_.pop_back();
    }

    void Triangulation::ensure_mesh_()
    {
        if (mesh_built_) return;

        mesh_vertices_.clear();
        mesh_halfedges_.clear();
        free_halfedges_.clear();
        last_face_he_ = 0;
        triangle_to_he_.clear();

        const size_t nTri = triangles.size();
        const size_t nHE = 3 * nTri;

        std::vector<int> orig_to_local(input_points.size(), -1);
        mesh_vertices_.reserve(nTri * 3);

        for (size_t t = 0; t < nTri; ++t)
        {
            const auto& tri = triangles[t];
            for (int i = 0; i < 3; ++i)
            {
                const size_t orig = tri[i];
                int& slot = orig_to_local[orig];
                if (slot < 0)
                {
                    slot = static_cast<int>(mesh_vertices_.size());
                    MeshVertex v;
                    v.p = input_points[orig];
                    v.orig_id = orig;
                    v.halfedge = -1;
                    mesh_vertices_.push_back(v);
                }
            }
        }

        mesh_halfedges_.resize(nHE);
        triangle_to_he_.resize(nTri);

        MeshHalfEdge* __restrict he_ptr = mesh_halfedges_.data();
        for (size_t t = 0; t < nTri; ++t)
        {
            const auto& tri = triangles[t];
            const int v0 = orig_to_local[tri[0]];
            const int v1 = orig_to_local[tri[1]];
            const int v2 = orig_to_local[tri[2]];

            const int base = static_cast<int>(3 * t);

            he_ptr[base + 0] = { v0, -1, base + 1, static_cast<int>(t), false };
            he_ptr[base + 1] = { v1, -1, base + 2, static_cast<int>(t), false };
            he_ptr[base + 2] = { v2, -1, base + 0, static_cast<int>(t), false };

            if (mesh_vertices_[v0].halfedge < 0) mesh_vertices_[v0].halfedge = base + 0;
            if (mesh_vertices_[v1].halfedge < 0) mesh_vertices_[v1].halfedge = base + 1;
            if (mesh_vertices_[v2].halfedge < 0) mesh_vertices_[v2].halfedge = base + 2;

            triangle_to_he_[t] = base;
        }

        const int nV = static_cast<int>(mesh_vertices_.size());

        std::vector<int> offset(nV + 1, 0);
        for (size_t i = 0; i < nHE; ++i)
            ++offset[he_ptr[i].origin + 1];
        for (int v = 0; v < nV; ++v)
            offset[v + 1] += offset[v];

        std::vector<int> adj(nHE);
        {
            std::vector<int> cursor(offset.begin(), offset.begin() + nV);
            for (size_t i = 0; i < nHE; ++i)
            {
                const int u = he_ptr[i].origin;
                adj[cursor[u]++] = static_cast<int>(i);
            }
        }

        for (size_t i = 0; i < nHE; ++i)
        {
            if (he_ptr[i].twin >= 0) continue;

            const int u = he_ptr[i].origin;
            const int v = he_ptr[he_ptr[i].next].origin;

            const int start = offset[v];
            const int end = offset[v + 1];

            for (int k = start; k < end; ++k)
            {
                const int h2 = adj[k];
                if (he_ptr[h2].origin != v) continue;
                const int u2 = he_ptr[he_ptr[h2].next].origin;
                if (u2 == u)
                {
                    he_ptr[i].twin = h2;
                    he_ptr[h2].twin = static_cast<int>(i);
                    break;
                }
            }
        }

        mesh_built_ = true;
    }

    void Triangulation::ensure_mesh_parallel_()
    {
        if (mesh_built_) return;

        unsigned nthreads = std::thread::hardware_concurrency();
        if (nthreads == 0) nthreads = 1;
        if (nthreads > 64) nthreads = 64;

        const size_t nTri = triangles.size();

        if (nthreads == 1 || nTri < 50000)
        {
            ensure_mesh_();
            return;
        }

        mesh_vertices_.clear();
        mesh_halfedges_.clear();
        free_halfedges_.clear();
        last_face_he_ = 0;
        triangle_to_he_.clear();

        const size_t nHE = 3 * nTri;

        std::vector<int> orig_to_local(input_points.size(), -1);
        mesh_vertices_.reserve(nTri * 3);

        for (size_t t = 0; t < nTri; ++t)
        {
            const auto& tri = triangles[t];
            for (int i = 0; i < 3; ++i)
            {
                const size_t orig = tri[i];
                int& slot = orig_to_local[orig];
                if (slot < 0)
                {
                    slot = static_cast<int>(mesh_vertices_.size());
                    MeshVertex v;
                    v.p = input_points[orig];
                    v.orig_id = orig;
                    v.halfedge = -1;
                    mesh_vertices_.push_back(v);
                }
            }
        }

        mesh_halfedges_.resize(nHE);
        triangle_to_he_.resize(nTri);

        const int nV = static_cast<int>(mesh_vertices_.size());

        std::unique_ptr<std::atomic<int>[]> v_he(new std::atomic<int>[nV]);
        for (int i = 0; i < nV; ++i)
            v_he[i].store(-1, std::memory_order_relaxed);

        {
            std::vector<std::thread> workers;
            workers.reserve(nthreads);

            auto fill = [&](size_t t0, size_t t1) {
                MeshHalfEdge* hp = mesh_halfedges_.data();
                for (size_t t = t0; t < t1; ++t)
                {
                    const auto& tri = triangles[t];
                    const int v0 = orig_to_local[tri[0]];
                    const int v1 = orig_to_local[tri[1]];
                    const int v2 = orig_to_local[tri[2]];
                    const int base = static_cast<int>(3 * t);

                    hp[base + 0] = { v0, -1, base + 1, static_cast<int>(t), false };
                    hp[base + 1] = { v1, -1, base + 2, static_cast<int>(t), false };
                    hp[base + 2] = { v2, -1, base + 0, static_cast<int>(t), false };

                    triangle_to_he_[t] = base;

                    int expected = -1;
                    v_he[v0].compare_exchange_strong(
                        expected, base + 0, std::memory_order_relaxed);
                    expected = -1;
                    v_he[v1].compare_exchange_strong(
                        expected, base + 1, std::memory_order_relaxed);
                    expected = -1;
                    v_he[v2].compare_exchange_strong(
                        expected, base + 2, std::memory_order_relaxed);
                }
                };

            const size_t chunk = (nTri + nthreads - 1) / nthreads;
            for (unsigned i = 0; i < nthreads; ++i)
            {
                const size_t t0 = i * chunk;
                const size_t t1 = std::min(nTri, t0 + chunk);
                if (t0 >= t1) break;
                workers.emplace_back(fill, t0, t1);
            }
            for (auto& th : workers) th.join();
        }

        for (int i = 0; i < nV; ++i)
            mesh_vertices_[i].halfedge = v_he[i].load(std::memory_order_relaxed);

        std::vector<int> offset(nV + 1, 0);
        {
            MeshHalfEdge* hp = mesh_halfedges_.data();
            for (size_t i = 0; i < nHE; ++i)
                ++offset[hp[i].origin + 1];
            for (int v = 0; v < nV; ++v)
                offset[v + 1] += offset[v];
        }

        std::vector<int> adj(nHE);
        {
            MeshHalfEdge* hp = mesh_halfedges_.data();
            std::vector<int> cursor(offset.begin(), offset.begin() + nV);
            for (size_t i = 0; i < nHE; ++i)
            {
                const int u = hp[i].origin;
                adj[cursor[u]++] = static_cast<int>(i);
            }
        }

        {
            std::vector<std::thread> workers;
            workers.reserve(nthreads);

            auto find_twins = [&](size_t i0, size_t i1) {
                MeshHalfEdge* hp = mesh_halfedges_.data();
                for (size_t i = i0; i < i1; ++i)
                {
                    if (hp[i].twin >= 0) continue;

                    const int u = hp[i].origin;
                    const int v = hp[hp[i].next].origin;

                    const int start = offset[v];
                    const int end = offset[v + 1];
                    for (int k = start; k < end; ++k)
                    {
                        const int h2 = adj[k];
                        if (hp[h2].origin != v) continue;
                        const int u2 = hp[hp[h2].next].origin;
                        if (u2 == u)
                        {
                            hp[i].twin = h2;
                            break;
                        }
                    }
                }
                };

            const size_t chunk = (nHE + nthreads - 1) / nthreads;
            for (unsigned i = 0; i < nthreads; ++i)
            {
                const size_t i0 = i * chunk;
                const size_t i1 = std::min(nHE, i0 + chunk);
                if (i0 >= i1) break;
                workers.emplace_back(find_twins, i0, i1);
            }
            for (auto& th : workers) th.join();
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
        mesh_halfedges_[h].triangle_idx = -1;
        mesh_halfedges_[h].is_constrained = false;
        free_halfedges_.push_back(h);
    }

    int Triangulation::locate_(const qdt_vec2d& p, int& hull_he)
    {
        hull_he = -1;
        if (mesh_halfedges_.empty()) return -1;

        int h = last_face_he_;
        if (h < 0 || h >= static_cast<int>(mesh_halfedges_.size()) ||
            mesh_halfedges_[h].origin < 0)
        {
            h = 0;
            while (h < static_cast<int>(mesh_halfedges_.size()) &&
                mesh_halfedges_[h].origin < 0)
                ++h;
            if (h >= static_cast<int>(mesh_halfedges_.size())) return -1;
        }

        // A walk that is not stuck visits every triangle at most once. A
        // triangulation that is not Delaunay (constrained edges that were
        // flipped into place) can make the plain visibility walk run in
        // circles, so after this many steps the point is looked up by
        // scanning instead.
        const size_t max_steps = triangle_to_he_.size() + 16;
        for (size_t step = 0; step < max_steps; ++step)
        {
            const int h0 = h;
            const int h1 = mesh_halfedges_[h0].next;
            const int h2 = mesh_halfedges_[h1].next;
            if (mesh_halfedges_[h2].next != h0) break;

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

        return locate_scan_(p, hull_he);
    }

    // Finds the triangle that contains p (boundary included) by testing all of
    // them. Returns a halfedge of that triangle. If p is outside the convex
    // hull, returns -1 and sets hull_he to a hull halfedge that p can see.
    int Triangulation::locate_scan_(const qdt_vec2d& p, int& hull_he)
    {
        hull_he = -1;
        for (size_t t = 0; t < triangle_to_he_.size(); ++t)
        {
            const int h0 = triangle_to_he_[t];
            const int h1 = mesh_halfedges_[h0].next;
            const int h2 = mesh_halfedges_[h1].next;
            const qdt_vec2d& A = mesh_vertices_[mesh_halfedges_[h0].origin].p;
            const qdt_vec2d& B = mesh_vertices_[mesh_halfedges_[h1].origin].p;
            const qdt_vec2d& C = mesh_vertices_[mesh_halfedges_[h2].origin].p;
            if (robust::orient2d(A, B, p) >= 0 && robust::orient2d(B, C, p) >= 0 &&
                robust::orient2d(C, A, p) >= 0)
            {
                last_face_he_ = h0;
                return h0;
            }
        }
        for (size_t h = 0; h < mesh_halfedges_.size(); ++h)
        {
            const MeshHalfEdge& e = mesh_halfedges_[h];
            if (e.origin < 0 || e.twin >= 0) continue;
            if (robust::orient2d(mesh_vertices_[e.origin].p,
                mesh_vertices_[mesh_halfedges_[e.next].origin].p, p) < 0)
            {
                hull_he = static_cast<int>(h);
                return -1;
            }
        }
        return -1;
    }

    void Triangulation::insert_by_rebuild_(const qdt_vec2d& p)
    {
        input_points.push_back(p);
        Triangulation fresh = delaunay_triangulate(
            input_points.data(), input_points.size());
        *this = std::move(fresh);
    }

    void Triangulation::insert(const qdt_vec2d& p)
    {
        ensure_mesh_();
        insert_after_mesh_(p);
    }

    void Triangulation::insert_parallel(const qdt_vec2d& p)
    {
        ensure_mesh_parallel_();
        insert_after_mesh_(p);
    }

    void Triangulation::insert_after_mesh_(const qdt_vec2d& p)
    {
        for (const auto& v : mesh_vertices_)
        {
            if (v.p.x == p.x && v.p.y == p.y)
            {
                input_points.push_back(p);
                rep.push_back(v.orig_id);
                return;
            }
        }

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
        const bool outside = (h_start < 0);

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

        std::vector<int> visible_hull;
        if (outside)
        {
            auto is_visible = [&](int e) {
                const qdt_vec2d& A = mesh_vertices_[mesh_halfedges_[e].origin].p;
                const qdt_vec2d& B = mesh_vertices_[he_dest(e)].p;
                return robust::orient2d(A, B, p) < 0;
                };
            auto next_hull = [&](int e) {
                int x = mesh_halfedges_[e].next;
                while (mesh_halfedges_[x].twin >= 0)
                    x = mesh_halfedges_[mesh_halfedges_[x].twin].next;
                return x;
                };
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

        // The cavity is only the triangle that contains p, plus the triangle on
        // the other side of the edge p lies on (if any), so that edge is split
        // and no degenerate triangle appears. Growing the cavity with
        // circumcircle tests alone would also swallow constrained edges inside
        // the grown region and delete them. Delaunay is restored afterwards by
        // edge flips (see the end of this function), which never flip a
        // constrained edge, so constraints survive.
        // If p is outside the hull nothing is removed: p is joined to the
        // visible hull edges and the flips take care of the rest.
        std::vector<int>  cavity_hes;
        std::vector<bool> visited_tri(triangles.size(), false);
        int split_a = -1, split_b = -1; // endpoints of the constrained edge p splits

        if (!outside)
        {
            auto add_to_cavity = [&](int h) {
                const int t = mesh_halfedges_[h].triangle_idx;
                if (t >= 0 && !visited_tri[t])
                {
                    visited_tri[t] = true;
                    cavity_hes.push_back(h);
                }
                };

            add_to_cavity(h_start);
            int e = h_start;
            do
            {
                const int e_next = mesh_halfedges_[e].next;
                if (robust::orient2d(mesh_vertices_[mesh_halfedges_[e].origin].p,
                    mesh_vertices_[mesh_halfedges_[e_next].origin].p, p) == 0)
                {
                    // p lies on this edge, so the triangle across it is split too.
                    const int tw = mesh_halfedges_[e].twin;
                    if (tw >= 0) add_to_cavity(tw);
                    if (mesh_halfedges_[e].is_constrained)
                    {
                        split_a = mesh_halfedges_[e].origin;
                        split_b = mesh_halfedges_[e_next].origin;
                    }
                }
                e = e_next;
            } while (e != h_start);
        }

        struct BoundaryInfo { int a, b, external_twin; bool is_constrained; };
        std::vector<BoundaryInfo> boundary_info;

        for (int h : cavity_hes)
        {
            int e = h;
            do
            {
                const int tw = mesh_halfedges_[e].twin;
                int tw_tri = (tw >= 0) ? mesh_halfedges_[tw].triangle_idx : -1;
                bool tw_in_cavity = (tw_tri >= 0 && visited_tri[tw_tri]);

                if (!tw_in_cavity)
                {
                    const int a = mesh_halfedges_[e].origin;
                    const int b = mesh_halfedges_[mesh_halfedges_[e].next].origin;
                    if (tw >= 0 ||
                        robust::orient2d(mesh_vertices_[a].p, mesh_vertices_[b].p, p) > 0)
                        boundary_info.push_back({ a, b, tw, mesh_halfedges_[e].is_constrained });
                }

                e = mesh_halfedges_[e].next;
            } while (e != h);
        }

        for (int e : visible_hull)
        {
            const int tri_idx = mesh_halfedges_[e].triangle_idx;
            if (tri_idx >= 0 && visited_tri[tri_idx])
                continue;
            const int e1 = mesh_halfedges_[e].next;
            const int a = mesh_halfedges_[e].origin;
            const int b = mesh_halfedges_[e1].origin;
            boundary_info.push_back({ b, a, e, mesh_halfedges_[e].is_constrained });
        }

        std::vector<int> cavity_tri_indices;
        for (size_t i = 0; i < visited_tri.size(); ++i)
        {
            if (visited_tri[i])
                cavity_tri_indices.push_back(static_cast<int>(i));
        }
        std::sort(cavity_tri_indices.rbegin(), cavity_tri_indices.rend());
        for (int idx : cavity_tri_indices)
        {
            remove_triangle_by_index_(idx);
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

            mesh_halfedges_[h0] = { bi.a, bi.external_twin, h1, -1, bi.is_constrained };
            mesh_halfedges_[h1] = { bi.b, -1, h2, -1, false };
            mesh_halfedges_[h2] = { new_v, -1, h0, -1, false };

            if (bi.external_twin >= 0)
            {
                mesh_halfedges_[bi.external_twin].twin = h0;
                mesh_halfedges_[bi.external_twin].is_constrained = bi.is_constrained;
            }

            size_t new_tri_idx = add_triangle_to_output_(
                mesh_vertices_[bi.a].orig_id,
                mesh_vertices_[bi.b].orig_id,
                new_id,
                h0
            );

            mesh_halfedges_[h0].triangle_idx = static_cast<int>(new_tri_idx);
            mesh_halfedges_[h1].triangle_idx = static_cast<int>(new_tri_idx);
            mesh_halfedges_[h2].triangle_idx = static_cast<int>(new_tri_idx);

            new_tris.push_back({ h0, h1, h2 });
        }

        if (new_tris.empty())
        {
            Triangulation fresh = delaunay_triangulate(
                input_points.data(), input_points.size());
            *this = std::move(fresh);
            return;
        }

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

        // The removed triangles' halfedges were freed and their slots reused
        // above, so the stored halfedge of every cavity vertex may be stale.
        // h0 starts at boundary_info[i].a and h1 at boundary_info[i].b.
        for (size_t i = 0; i < new_tris.size(); ++i)
        {
            mesh_vertices_[boundary_info[i].a].halfedge = new_tris[i].h0;
            mesh_vertices_[boundary_info[i].b].halfedge = new_tris[i].h1;
        }

        // p split a constrained edge (split_a, split_b): both halves stay
        // constrained. new_tris[i] was built from boundary_info[i]; its h1
        // runs b -> p and its h2 runs p -> a.
        if (split_a >= 0)
        {
            auto mark_constrained = [&](int h) {
                mesh_halfedges_[h].is_constrained = true;
                const int tw = mesh_halfedges_[h].twin;
                if (tw >= 0) mesh_halfedges_[tw].is_constrained = true;
                };
            for (size_t i = 0; i < new_tris.size(); ++i)
            {
                if (boundary_info[i].b == split_a || boundary_info[i].b == split_b)
                    mark_constrained(new_tris[i].h1);
                if (boundary_info[i].a == split_a || boundary_info[i].a == split_b)
                    mark_constrained(new_tris[i].h2);
            }
        }

        // Restore the (constrained) Delaunay property around p with edge
        // flips. h0 of every new triangle is the edge opposite p; flipping it
        // puts a diagonal through p, and the two edges that then face p are
        // checked in turn. Constrained edges are never flipped.
        std::vector<int> pending;
        for (const auto& nt : new_tris)
            pending.push_back(nt.h0);
        while (!pending.empty())
        {
            const int h = pending.back();
            pending.pop_back();

            const int tw = mesh_halfedges_[h].twin;
            if (tw < 0 || mesh_halfedges_[h].is_constrained)
                continue;

            const int h_next = mesh_halfedges_[h].next;
            const int h_prev = mesh_halfedges_[h_next].next;
            if (mesh_halfedges_[h_prev].origin != new_v)
                continue; // h's triangle no longer has p as its apex

            const int tw_next = mesh_halfedges_[tw].next;
            const int tw_prev = mesh_halfedges_[tw_next].next;

            const int a = mesh_halfedges_[h].origin;
            const int b = mesh_halfedges_[h_next].origin;
            const int d = mesh_halfedges_[tw_prev].origin;

            if (robust::incircle(mesh_vertices_[a].p, mesh_vertices_[b].p,
                mesh_vertices_[new_v].p, mesh_vertices_[d].p) <= 0)
                continue;

            if (flip_edge_(h))
            {
                pending.push_back(tw_next); // a -> d
                pending.push_back(tw_prev); // d -> b
            }
        }
    }

    int Triangulation::orig_to_local_vtx_(size_t orig_id) const
    {
        for (size_t i = 0; i < mesh_vertices_.size(); ++i)
        {
            if (mesh_vertices_[i].orig_id == orig_id)
                return static_cast<int>(i);
        }
        return -1;
    }

    // Calls f(h) for every halfedge h that starts at vertex x, until f returns
    // true. Works for interior and for convex hull vertices (whose fan is open
    // on the outer side), independent of which halfedge is stored for x.
    template <class F>
    void Triangulation::for_each_out_(int x, F&& f) const
    {
        const int h0 = mesh_vertices_[x].halfedge;
        if (h0 < 0) return;
        const size_t limit = mesh_halfedges_.size();
        size_t steps = 0;

        // Counter-clockwise from h0: the next halfedge out of x is the twin of
        // the previous halfedge of the same triangle.
        int h = h0;
        for (;;)
        {
            if (f(h)) return;
            if (++steps > limit) return;
            const int prev = mesh_halfedges_[mesh_halfedges_[h].next].next;
            const int tw = mesh_halfedges_[prev].twin;
            if (tw < 0) break;      // reached the hull boundary
            if (tw == h0) return;   // went all the way around
            h = tw;
        }

        // Clockwise from h0 (only reached for a hull vertex).
        h = h0;
        for (;;)
        {
            const int tw = mesh_halfedges_[h].twin;
            if (tw < 0) return;
            h = mesh_halfedges_[tw].next;
            if (h == h0) return;
            if (f(h)) return;
            if (++steps > limit) return;
        }
    }

    // The halfedge that goes from u to v, or -1. A convex hull edge has a
    // halfedge in one direction only, so callers that accept either direction
    // look for (u, v) and then for (v, u).
    int Triangulation::find_halfedge_(int u_local, int v_local) const
    {
        int found = -1;
        for_each_out_(u_local, [&](int h) -> bool {
            if (mesh_halfedges_[mesh_halfedges_[h].next].origin == v_local)
            {
                found = h;
                return true;
            }
            return false;
            });
        return found;
    }

    bool Triangulation::flip_edge_(int h)
    {
        int tw = mesh_halfedges_[h].twin;
        if (tw < 0) return false;
        if (mesh_halfedges_[h].is_constrained) return false;

        int h_next = mesh_halfedges_[h].next;
        int h_prev = mesh_halfedges_[h_next].next;

        int tw_next = mesh_halfedges_[tw].next;
        int tw_prev = mesh_halfedges_[tw_next].next;

        // Before the flip (both triangles counter-clockwise):
        //   triangle 1:  h : a->b,  h_next : b->c,  h_prev : c->a
        //   triangle 2:  tw: b->a,  tw_next: a->d,  tw_prev: d->b
        // so c is the apex of triangle 1 and d is the apex of triangle 2.
        // (tw originates at b and tw_next originates at a, so neither can be
        // used for c or d.)
        const int a = mesh_halfedges_[h].origin;
        const int b = mesh_halfedges_[h_next].origin;
        const int c = mesh_halfedges_[h_prev].origin;
        const int d = mesh_halfedges_[tw_prev].origin;

        // After the flip the diagonal is c-d and the triangles are (d,c,a) and
        // (c,d,b). The flip is only valid when both stay strictly
        // counter-clockwise, i.e. the quadrilateral a,d,b,c is strictly convex.
        if (robust::orient2d(mesh_vertices_[d].p, mesh_vertices_[c].p, mesh_vertices_[a].p) <= 0 ||
            robust::orient2d(mesh_vertices_[c].p, mesh_vertices_[d].p, mesh_vertices_[b].p) <= 0)
        {
            return false;
        }

        // triangle 1 becomes  h : d->c,  h_prev : c->a,  tw_next : a->d
        mesh_halfedges_[h].origin = d;
        mesh_halfedges_[h].next = h_prev;
        mesh_halfedges_[h_prev].next = tw_next;
        mesh_halfedges_[tw_next].next = h;

        // triangle 2 becomes  tw : c->d,  tw_prev : d->b,  h_next : b->c
        mesh_halfedges_[tw].origin = c;
        mesh_halfedges_[tw].next = tw_prev;
        mesh_halfedges_[tw_prev].next = h_next;
        mesh_halfedges_[h_next].next = tw;

        const int tri1 = mesh_halfedges_[h].triangle_idx;
        const int tri2 = mesh_halfedges_[tw].triangle_idx;

        // h and h_prev stay in triangle 1, tw and tw_prev stay in triangle 2;
        // tw_next moves to triangle 1 and h_next moves to triangle 2.
        mesh_halfedges_[tw_next].triangle_idx = tri1;
        mesh_halfedges_[h_next].triangle_idx = tri2;

        triangles[tri1] = { mesh_vertices_[d].orig_id, mesh_vertices_[c].orig_id, mesh_vertices_[a].orig_id };
        triangles[tri2] = { mesh_vertices_[c].orig_id, mesh_vertices_[d].orig_id, mesh_vertices_[b].orig_id };

        // remove_triangle_by_index_ walks the cycle starting at
        // triangle_to_he_[i], so the representative must lie in triangle i.
        triangle_to_he_[tri1] = h;
        triangle_to_he_[tri2] = tw;

        // A vertex's stored halfedge must start at that vertex.
        mesh_vertices_[a].halfedge = tw_next;
        mesh_vertices_[b].halfedge = h_next;
        mesh_vertices_[c].halfedge = h_prev;
        mesh_vertices_[d].halfedge = h;

        return true;
    }

    bool Triangulation::segment_intersects_(const qdt_vec2d& a, const qdt_vec2d& b, const qdt_vec2d& c, const qdt_vec2d& d) const
    {
        double o1 = robust::orient2d(a, b, c);
        double o2 = robust::orient2d(a, b, d);
        double o3 = robust::orient2d(c, d, a);
        double o4 = robust::orient2d(c, d, b);

        return ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
            ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0));
    }

    void Triangulation::mark_constrained_(int h)
    {
        mesh_halfedges_[h].is_constrained = true;
        const int tw = mesh_halfedges_[h].twin;
        if (tw >= 0) mesh_halfedges_[tw].is_constrained = true;
    }

    // Do the vectors a->b and a->c point the same way? Only meaningful when
    // they are known to be collinear. Compares the signs of the coordinate
    // differences, so it cannot underflow.
    static bool same_direction_(const qdt_vec2d& a, const qdt_vec2d& b, const qdt_vec2d& c)
    {
        auto sg = [](double d) { return d > 0 ? 1 : (d < 0 ? -1 : 0); };
        const int sx = sg(b.x - a.x) * sg(c.x - a.x);
        const int sy = sg(b.y - a.y) * sg(c.y - a.y);
        return sx >= 0 && sy >= 0 && (sx > 0 || sy > 0);
    }

    // Walks from u to v and lists, in order, the edges the open segment uv
    // crosses. Each is returned as the halfedge that runs from the vertex on
    // the right of u->v to the vertex on its left. Nothing is modified.
    //
    // Returns false (and the list is meaningless) when the segment cannot
    // become an edge as it is: it crosses a constrained edge, or it passes
    // exactly through another vertex.
    bool Triangulation::trace_segment_(int u, int v, std::vector<int>& crossed) const
    {
        auto he_next = [&](int h) { return mesh_halfedges_[h].next; };
        auto org = [&](int h) { return mesh_halfedges_[h].origin; };
        auto pt = [&](int x) -> const qdt_vec2d& { return mesh_vertices_[x].p; };
        const qdt_vec2d& U = pt(u);
        const qdt_vec2d& V = pt(v);

        // Among the triangles around u, find the one the segment leaves u
        // through: v is strictly between its edges u->y and u->z. If v lies
        // exactly on one of those edges, that edge's other end is a vertex on
        // the segment.
        int start = -1; // y -> z, the edge opposite u in that triangle
        bool blocked = false;
        for_each_out_(u, [&](int h) -> bool {
            const int hn = he_next(h);
            const int y = org(hn);
            const int z = org(he_next(hn));
            const double oy = robust::orient2d(U, pt(y), V);
            const double oz = robust::orient2d(U, pt(z), V);
            if (oy > 0 && oz < 0)
            {
                start = hn;
                return true;
            }
            if ((oy == 0 && same_direction_(U, pt(y), V)) ||
                (oz == 0 && same_direction_(U, pt(z), V)))
            {
                blocked = true;
                return true;
            }
            return false;
            });
        if (blocked || start < 0) return false;

        // In the first triangle y lies on the right of u->v and z on its left.
        // Every crossed edge keeps that orientation (right vertex -> left vertex).
        int c = start;
        for (size_t guard = 0; guard <= triangles.size(); ++guard)
        {
            if (mesh_halfedges_[c].is_constrained) return false;
            crossed.push_back(c);

            const int tw = mesh_halfedges_[c].twin; // z -> y in the next triangle
            if (tw < 0) return false;               // would leave the convex hull
            const int a = he_next(tw);              // y -> w
            const int b = he_next(a);               // w -> z
            const int w = org(b);
            if (w == v) return true;

            const double o = robust::orient2d(U, V, pt(w));
            if (o == 0) return false;               // w lies on the segment
            c = (o > 0) ? a : b;                    // w on the left: leave via y->w
        }
        return false;
    }

    // Lawson's flip algorithm: every edge on the stack that is not
    // constrained and not locally Delaunay is flipped, and the four edges
    // around the flipped one are examined in turn.
    void Triangulation::legalize_(std::vector<int>& stack)
    {
        while (!stack.empty())
        {
            const int h = stack.back();
            stack.pop_back();

            if (mesh_halfedges_[h].origin < 0) continue;
            const int tw = mesh_halfedges_[h].twin;
            if (tw < 0 || mesh_halfedges_[h].is_constrained) continue;

            const int h_next = mesh_halfedges_[h].next;
            const int h_prev = mesh_halfedges_[h_next].next;
            const int tw_next = mesh_halfedges_[tw].next;
            const int tw_prev = mesh_halfedges_[tw_next].next;

            const int a = mesh_halfedges_[h].origin;
            const int b = mesh_halfedges_[h_next].origin;
            const int c = mesh_halfedges_[h_prev].origin;
            const int d = mesh_halfedges_[tw_prev].origin;

            if (robust::incircle(mesh_vertices_[a].p, mesh_vertices_[b].p,
                mesh_vertices_[c].p, mesh_vertices_[d].p) <= 0)
                continue;
            if (!flip_edge_(h)) continue;

            stack.push_back(h_prev);
            stack.push_back(h_next);
            stack.push_back(tw_next);
            stack.push_back(tw_prev);
        }
    }

    // Inserts the constraint between two mesh vertices (local indices).
    //  - If the edge exists it just becomes constrained.
    //  - Otherwise the edges the segment crosses are flipped away (Sloan 1993:
    //    a crossed edge whose quadrilateral is strictly convex is flipped; if
    //    the new diagonal still crosses the segment it goes back in the queue),
    //    which always ends with the edge u-v in place, and then the
    //    constrained Delaunay property is restored with Lawson flips.
    //  - If the segment crosses a constrained edge or passes through another
    //    vertex it cannot be an edge: nothing is changed and false is returned.
    bool Triangulation::insert_constraint_local_(int u, int v)
    {
        int e = find_halfedge_(u, v);
        if (e < 0) e = find_halfedge_(v, u);
        if (e >= 0)
        {
            mark_constrained_(e);
            return true;
        }

        std::vector<int> crossed;
        if (!trace_segment_(u, v, crossed)) return false;

        const qdt_vec2d& U = mesh_vertices_[u].p;
        const qdt_vec2d& V = mesh_vertices_[v].p;

        std::deque<int> queue(crossed.begin(), crossed.end());
        std::vector<int> new_edges;
        size_t stalled = 0;
        while (!queue.empty())
        {
            const int h = queue.front();
            queue.pop_front();
            if (flip_edge_(h))
            {
                stalled = 0;
                const int c = mesh_halfedges_[h].origin;
                const int d = mesh_halfedges_[mesh_halfedges_[h].next].origin;
                if (segment_intersects_(U, V, mesh_vertices_[c].p, mesh_vertices_[d].p))
                    queue.push_back(h);
                else
                    new_edges.push_back(h);
            }
            else
            {
                queue.push_back(h);
                // A whole round without a single flip. With no vertex on the
                // segment this cannot happen; stop instead of looping.
                if (++stalled >= queue.size()) break;
            }
        }

        bool ok = false;
        if (queue.empty())
        {
            e = -1;
            for (int h : new_edges)
            {
                const int a = mesh_halfedges_[h].origin;
                const int b = mesh_halfedges_[mesh_halfedges_[h].next].origin;
                if ((a == u && b == v) || (a == v && b == u))
                {
                    e = h;
                    break;
                }
            }
            if (e < 0) e = find_halfedge_(u, v);
            if (e < 0) e = find_halfedge_(v, u);
            if (e >= 0)
            {
                mark_constrained_(e);
                ok = true;
            }
        }

        // Restore the Delaunay property around everything that changed. The
        // constrained edge u-v is already marked, so it is never flipped.
        std::vector<int> stack;
        auto push_around = [&](int h) {
            if (mesh_halfedges_[h].origin < 0) return;
            const int n = mesh_halfedges_[h].next;
            stack.push_back(h);
            stack.push_back(n);
            stack.push_back(mesh_halfedges_[n].next);
            const int tw = mesh_halfedges_[h].twin;
            if (tw >= 0)
            {
                const int tn = mesh_halfedges_[tw].next;
                stack.push_back(tn);
                stack.push_back(mesh_halfedges_[tn].next);
            }
            };
        for (int h : new_edges) push_around(h);
        for (int h : queue) push_around(h);
        legalize_(stack);
        return ok;
    }

    bool Triangulation::insert_constraint(size_t orig_v0, size_t orig_v1)
    {
        // Indices come from the caller: reject anything that is not a point of
        // this context (this also covers a context that was never filled).
        if (orig_v0 >= rep.size() || orig_v1 >= rep.size()) return false;
        orig_v0 = rep[orig_v0];
        orig_v1 = rep[orig_v1];
        if (orig_v0 == orig_v1) return false;

        ensure_mesh_();
        const int u = orig_to_local_vtx_(orig_v0);
        const int v = orig_to_local_vtx_(orig_v1);
        if (u < 0 || v < 0) return false;
        return insert_constraint_local_(u, v);
    }

    size_t Triangulation::insert_constraints_batch(const qdt_edge* edges, size_t num_edges)
    {
        if (!edges || num_edges == 0) return 0;
        ensure_mesh_();

        // Constraints never add vertices, so one index table serves the whole
        // batch instead of one linear search per endpoint.
        std::vector<int> to_local(input_points.size(), -1);
        for (size_t i = 0; i < mesh_vertices_.size(); ++i)
        {
            if (mesh_vertices_[i].orig_id < to_local.size())
                to_local[mesh_vertices_[i].orig_id] = static_cast<int>(i);
        }

        size_t count = 0;
        for (size_t i = 0; i < num_edges; ++i)
        {
            size_t a = edges[i].v0;
            size_t b = edges[i].v1;
            if (a >= rep.size() || b >= rep.size()) continue;
            a = rep[a];
            b = rep[b];
            if (a == b || a >= to_local.size() || b >= to_local.size()) continue;
            const int u = to_local[a];
            const int v = to_local[b];
            if (u < 0 || v < 0) continue;
            if (insert_constraint_local_(u, v))
                ++count;
        }
        return count;
    }

    Triangulation delaunay_triangulate(Point* pts, size_t sz)
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
            detail::Builder b(V);
            b.run();
            b.extract(res.triangles);
        }

        res.input_points.resize(sz);
        for (size_t i = 0; i < sz; ++i)
            res.input_points[i] = pts[i];

        return res;
    }
}