/*
 * qdt_verify - correctness test for libqdt that can write its result into the
 * README.
 *
 * The program uses only the public C API (qdt.h). Every triangulation it
 * produces is checked with exact integer arithmetic (the test inputs have
 * integer coordinates with |c| <= 2^28, which doubles represent exactly), so
 * the verdict does not depend on the predicates inside the library.
 *
 * Usage:
 *   qdt_verify [options]
 *     --readme FILE        update FILE: the section between the markers
 *                          <!-- QDT-VERIFY-BEGIN --> and <!-- QDT-VERIFY-END -->
 *                          is replaced; the first time (no markers yet) the
 *                          section is inserted before "## Current limitations".
 *                          Only done when every scenario passed.
 *     --output FILE        write the Markdown section to FILE (UTF-8, LF)
 *     --note TEXT          add a line to the "Test environment" list (repeatable)
 *     --quick              fewer cases, no inputs with hundreds of thousands of
 *                          points (a smoke run)
 *     --iterations N       base number of random cases per scenario (default 3000)
 *     --seed N             offset of all random seeds (default 0)
 *     --scenario LIST      comma separated scenario names, or "all"
 *     --allow-partial      allow --readme together with --quick or --scenario
 *     --list               list the scenarios and exit
 *   Without --readme and --output the section is printed to stdout. Progress
 *   goes to stderr. Exit status: 0 all passed, 1 a check failed, 2 usage or
 *   file error.
 *
 * Build (see README, "Reproducing the benchmark"):
 *   g++ -std=c++14 -O2 -pthread -Ilibqdt test/qdt_verify.cpp libqdt.a -lm -o qdt_verify
 *
 * Define QDT_VERIFY_PORTABLE to use the built-in 160-bit integer type instead
 * of __int128 (this is what MSVC always uses; it is slower).
 */

#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS 1
#endif
#if defined(_WIN32)
#if !defined(NOMINMAX)
#define NOMINMAX 1
#endif
#if !defined(WIN32_LEAN_AND_MEAN)
#define WIN32_LEAN_AND_MEAN 1
#endif
#endif

#include "qdt.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#if defined(__linux__)
#include <sys/utsname.h>
#define QDT_VERIFY_LINUX 1
#endif
#if defined(_WIN32)
#include <windows.h>
#endif

 /* ------------------------------------------------------------------------ */
 /* Exact arithmetic                                                         */
 /* ------------------------------------------------------------------------ */

typedef long long ll;

#if defined(__SIZEOF_INT128__) && !defined(QDT_VERIFY_PORTABLE)

__extension__ typedef __int128 Wide;

static int sgn(const Wide& v)
{
    return v > 0 ? 1 : (v < 0 ? -1 : 0);
}

static const char* arith_name()
{
    return "`__int128`";
}

#else

/* 160-bit two's complement integer. Enough for the in-circle determinant of
   coordinates up to 2^28 (about 2^121). */
struct Wide
{
    uint32_t w[5];

    Wide()
    {
        for (int i = 0; i < 5; ++i)
        {
            w[i] = 0;
        }
    }

    Wide(long long v)
    {
        const uint64_t u = static_cast<uint64_t>(v);
        const uint32_t ext = v < 0 ? 0xFFFFFFFFu : 0u;
        w[0] = static_cast<uint32_t>(u);
        w[1] = static_cast<uint32_t>(u >> 32);
        w[2] = ext;
        w[3] = ext;
        w[4] = ext;
    }

    Wide& operator+=(const Wide& o)
    {
        uint64_t carry = 0;
        for (int i = 0; i < 5; ++i)
        {
            carry += static_cast<uint64_t>(w[i]) + o.w[i];
            w[i] = static_cast<uint32_t>(carry);
            carry >>= 32;
        }
        return *this;
    }
};

static Wide operator-(const Wide& a)
{
    Wide     r;
    uint64_t carry = 1;
    for (int i = 0; i < 5; ++i)
    {
        carry += static_cast<uint32_t>(~a.w[i]);
        r.w[i] = static_cast<uint32_t>(carry);
        carry >>= 32;
    }
    return r;
}

static Wide operator+(Wide a, const Wide& b)
{
    a += b;
    return a;
}

static Wide operator-(const Wide& a, const Wide& b)
{
    return a + (-b);
}

/* Truncated schoolbook product; exact modulo 2^160, which is enough here. */
static Wide operator*(const Wide& a, const Wide& b)
{
    Wide r;
    for (int i = 0; i < 5; ++i)
    {
        uint64_t carry = 0;
        for (int j = 0; i + j < 5; ++j)
        {
            const uint64_t cur = static_cast<uint64_t>(r.w[i + j]) + static_cast<uint64_t>(a.w[i]) * b.w[j] + carry;
            r.w[i + j] = static_cast<uint32_t>(cur);
            carry = cur >> 32;
        }
    }
    return r;
}

static int sgn(const Wide& a)
{
    if (a.w[4] & 0x80000000u)
    {
        return -1;
    }
    for (int i = 0; i < 5; ++i)
    {
        if (a.w[i])
        {
            return 1;
        }
    }
    return 0;
}

static const char* arith_name()
{
    return "built-in 160-bit integers";
}

#endif

struct P
{
    ll x, y;
};

typedef std::array<size_t, 3>        Tri;
typedef std::vector<Tri>             TriList;
typedef std::pair<size_t, size_t>    Edge;
typedef std::set<Edge>               ConSet;

static Wide orient(const P& a, const P& b, const P& c)
{
    return Wide(b.x - a.x) * Wide(c.y - a.y) - Wide(b.y - a.y) * Wide(c.x - a.x);
}

/* > 0 if d lies strictly inside the circle through a, b, c (counter-clockwise). */
static Wide incircle(const P& a, const P& b, const P& c, const P& d)
{
    const Wide adx(a.x - d.x), ady(a.y - d.y), bdx(b.x - d.x), bdy(b.y - d.y), cdx(c.x - d.x), cdy(c.y - d.y);
    const Wide al = adx * adx + ady * ady;
    const Wide bl = bdx * bdx + bdy * bdy;
    const Wide cl = cdx * cdx + cdy * cdy;
    return al * (bdx * cdy - cdx * bdy) + bl * (cdx * ady - adx * cdy) + cl * (adx * bdy - bdx * ady);
}

/* ------------------------------------------------------------------------ */
/* Small helpers                                                            */
/* ------------------------------------------------------------------------ */

typedef std::chrono::steady_clock Clock;

static double now_ms()
{
    return std::chrono::duration<double, std::milli>(Clock::now().time_since_epoch()).count();
}

static std::string fmt(const char* f, ...)
{
    char    buf[1024];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(buf, sizeof(buf), f, ap);
    va_end(ap);
    return std::string(buf);
}

static std::string with_commas(unsigned long long v)
{
    std::string s = std::to_string(v);
    for (int i = static_cast<int>(s.size()) - 3; i > 0; i -= 3)
    {
        s.insert(static_cast<size_t>(i), ",");
    }
    return s;
}

static void trim_str(std::string& s)
{
    size_t b = 0;
    while (b < s.size() && (s[b] == ' ' || s[b] == '\t' || s[b] == '\r' || s[b] == '\n'))
    {
        ++b;
    }
    size_t e = s.size();
    while (e > b && (s[e - 1] == ' ' || s[e - 1] == '\t' || s[e - 1] == '\r' || s[e - 1] == '\n'))
    {
        --e;
    }
    s = s.substr(b, e - b);
}

struct Config
{
    bool     quick = false;
    bool     partial = false; /* quick run or scenario subset */
    uint64_t iters = 3000;
    uint64_t seed = 0;
};

static uint64_t seed_of(const Config& c, uint64_t base)
{
    return base + c.seed * 1000003ULL;
}

/* Result of one scenario. pass() and fail() may be called from several threads. */
struct Scenario
{
    std::string              name;
    std::string              title;
    std::string              section;
    unsigned long long       cases = 0;
    unsigned long long       failures = 0;
    double                   ms = 0.0;
    bool                     skipped = false;
    std::string              skip_reason;
    std::string              note;
    std::vector<std::string> samples;
    std::mutex               mu;

    void pass()
    {
        std::lock_guard<std::mutex> g(mu);
        ++cases;
    }

    void fail(uint64_t seed, const std::string& msg)
    {
        std::lock_guard<std::mutex> g(mu);
        ++cases;
        ++failures;
        if (samples.size() < 5)
        {
            samples.push_back(fmt("seed %llu: ", static_cast<unsigned long long>(seed)) + msg);
            std::fprintf(stderr, "FAIL [%s] seed=%llu: %s\n", name.c_str(), static_cast<unsigned long long>(seed), msg.c_str());
        }
    }

    void expect(bool good, const char* what)
    {
        if (good)
        {
            pass();
        }
        else
        {
            fail(0, what);
        }
    }

    void finish(uint64_t seed, const std::string& msg)
    {
        if (msg.empty())
        {
            pass();
        }
        else
        {
            fail(seed, msg);
        }
    }
};

typedef void (*ExecFn)(qdt_t*, const qdt_vec2d*, size_t);
typedef size_t(*InsFn)(qdt_t*, qdt_vec2d);

static ExecFn exec_fn(bool par)
{
    return par ? qdt_execute_parallel : qdt_execute;
}

static InsFn ins_fn(bool par)
{
    return par ? qdt_insert_point_parallel : qdt_insert_point;
}

static std::vector<qdt_vec2d> to_d(const std::vector<P>& p, double scale = 1.0)
{
    std::vector<qdt_vec2d> r(p.size());
    for (size_t i = 0; i < p.size(); ++i)
    {
        r[i].x = static_cast<double>(p[i].x) * scale;
        r[i].y = static_cast<double>(p[i].y) * scale;
    }
    return r;
}

static std::vector<P> gen(std::mt19937_64& rng, size_t n, ll range)
{
    std::vector<P> p(n);
    for (size_t i = 0; i < n; ++i)
    {
        p[i].x = static_cast<ll>(rng() % static_cast<uint64_t>(range));
        p[i].y = static_cast<ll>(rng() % static_cast<uint64_t>(range));
    }
    return p;
}

/* Uniform in [-bound, bound]. */
static ll signed_rand(std::mt19937_64& rng, ll bound)
{
    return static_cast<ll>(rng() % static_cast<uint64_t>(2 * bound + 1)) - bound;
}

static TriList get_tris(const qdt_t* q, std::string* err)
{
    const size_t m = qdt_triangles_size(q);
    TriList      out(m);
    for (size_t i = 0; i < m; ++i)
    {
        if (!qdt_get_triangle(q, out[i].data(), i))
        {
            *err = "qdt_get_triangle failed for a valid id";
            return out;
        }
    }
    size_t sentinel[3] = { 7, 8, 9 };
    if (qdt_get_triangle(q, sentinel, m) != 0 || sentinel[0] != 7 || sentinel[1] != 8 || sentinel[2] != 9)
    {
        *err = "qdt_get_triangle(out of range) misbehaves";
    }
    return out;
}

/* Rotates a triangle so that its smallest index comes first (orientation kept). */
static Tri norm_tri(Tri t)
{
    int k = 0;
    if (t[1] < t[k]) k = 1;
    if (t[2] < t[k]) k = 2;
    const Tri r = { { t[k], t[(k + 1) % 3], t[(k + 2) % 3] } };
    return r;
}

/* Normalised and sorted: two triangulations are equal iff their lists are. */
static TriList tri_sorted(const TriList& v)
{
    TriList r(v.size());
    for (size_t i = 0; i < v.size(); ++i)
    {
        r[i] = norm_tri(v[i]);
    }
    std::sort(r.begin(), r.end());
    return r;
}

static bool coord_less(const P& a, const P& b)
{
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}

static bool coord_eq(const P& a, const P& b)
{
    return a.x == b.x && a.y == b.y;
}

/* Convex hull without collinear points, counter-clockwise. */
static std::vector<P> hull(std::vector<P> p)
{
    std::sort(p.begin(), p.end(), coord_less);
    p.erase(std::unique(p.begin(), p.end(), coord_eq), p.end());
    if (p.size() < 3)
    {
        return p;
    }
    std::vector<P> h(2 * p.size());
    size_t         k = 0;
    for (size_t i = 0; i < p.size(); ++i)
    {
        while (k >= 2 && sgn(orient(h[k - 2], h[k - 1], p[i])) <= 0) --k;
        h[k++] = p[i];
    }
    for (size_t i = p.size() - 1, t = k + 1; i > 0; --i)
    {
        while (k >= t && sgn(orient(h[k - 2], h[k - 1], p[i - 1])) <= 0) --k;
        h[k++] = p[i - 1];
    }
    h.resize(k - 1);
    return h;
}

/* True if the distinct points are at least three and not all collinear. */
static bool has_tris(const std::vector<P>& pts)
{
    std::vector<P> u = pts;
    std::sort(u.begin(), u.end(), coord_less);
    u.erase(std::unique(u.begin(), u.end(), coord_eq), u.end());
    if (u.size() < 3)
    {
        return false;
    }
    for (size_t i = 2; i < u.size(); ++i)
    {
        if (sgn(orient(u[0], u[1], u[i])) != 0)
        {
            return true;
        }
    }
    return false;
}

static bool on_open_segment(const P& a, const P& b, const P& p)
{
    if (sgn(orient(a, b, p)) != 0) return false;
    if (coord_eq(p, a) || coord_eq(p, b)) return false;
    return std::min(a.x, b.x) <= p.x && p.x <= std::max(a.x, b.x) && std::min(a.y, b.y) <= p.y &&
        p.y <= std::max(a.y, b.y);
}

static bool proper_cross(const P& a, const P& b, const P& c, const P& d)
{
    const int o1 = sgn(orient(a, b, c)), o2 = sgn(orient(a, b, d));
    const int o3 = sgn(orient(c, d, a)), o4 = sgn(orient(c, d, b));
    return o1 * o2 < 0 && o3 * o4 < 0;
}

/* ------------------------------------------------------------------------ */
/* The verifier                                                             */
/* ------------------------------------------------------------------------ */

struct DirEdge
{
    size_t a, b, c; /* directed edge a -> b of a triangle whose third vertex is c */
};

static bool de_less(const DirEdge& l, const DirEdge& r)
{
    return l.a < r.a || (l.a == r.a && l.b < r.b);
}

/* Finds the directed edge a -> b in the sorted list, or NULL. */
static const DirEdge* find_de(const std::vector<DirEdge>& es, size_t a, size_t b)
{
    DirEdge key = { a, b, 0 };
    auto    it = std::lower_bound(es.begin(), es.end(), key, de_less);
    if (it != es.end() && it->a == a && it->b == b)
    {
        return &*it;
    }
    return NULL;
}

/*
 * Checks the state of the context against the points it should hold.
 *   pts    all points of the context in index order (duplicates included)
 *   cons   constraints as (min, max) pairs of lowest-duplicate indices
 *   global also check the empty-circle property against every point
 *          (quadratic; for small inputs only)
 *   scale  the library holds pts[i] * scale
 * Returns "" if the context holds a valid (constrained) Delaunay triangulation
 * of exactly these points, otherwise a description of the first problem.
 */
static std::string verify(const qdt_t* q, const std::vector<P>& pts, const ConSet& cons, bool global, double scale = 1.0)
{
    std::string err;
    TriList     tris = get_tris(q, &err);
    if (!err.empty())
    {
        return err;
    }
    const size_t n = pts.size();

    if (qdt_points_size(q) != n)
    {
        return "qdt_points_size does not match the number of points";
    }
    if (qdt_point_is_used(q, n) != 0)
    {
        return "qdt_point_is_used(out of range) != 0";
    }
    for (size_t i = 0; i < n; ++i)
    {
        qdt_vec2d g = { 0.0, 0.0 };
        if (!qdt_get_point(q, &g, i) || g.x != static_cast<double>(pts[i].x) * scale ||
            g.y != static_cast<double>(pts[i].y) * scale)
        {
            return "qdt_get_point does not return the stored coordinates";
        }
    }
    {
        qdt_vec2d g = { -1.0, -1.0 };
        if (qdt_get_point(q, &g, n) != 0 || g.x != -1.0 || g.y != -1.0)
        {
            return "qdt_get_point(out of range) misbehaves";
        }
    }

    /* rep[i]: lowest index with the same coordinates as point i */
    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), static_cast<size_t>(0));
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (!coord_eq(pts[a], pts[b])) return coord_less(pts[a], pts[b]);
        return a < b;
        });
    std::vector<size_t> rep(n);
    for (size_t i = 0; i < n; ++i)
    {
        rep[order[i]] = (i > 0 && coord_eq(pts[order[i]], pts[order[i - 1]])) ? rep[order[i - 1]] : order[i];
    }
    std::vector<size_t> uniq;
    std::vector<P>      up;
    for (size_t i = 0; i < n; ++i)
    {
        if (rep[i] == i)
        {
            uniq.push_back(i);
            up.push_back(pts[i]);
        }
    }

    bool degenerate = uniq.size() < 3;
    if (!degenerate)
    {
        degenerate = true;
        for (size_t i = 2; i < uniq.size(); ++i)
        {
            if (sgn(orient(pts[uniq[0]], pts[uniq[1]], pts[uniq[i]])) != 0)
            {
                degenerate = false;
                break;
            }
        }
    }
    if (degenerate)
    {
        for (size_t i = 0; i < n; ++i)
        {
            if (qdt_point_is_used(q, i) != 0)
            {
                return "qdt_point_is_used != 0 although there are no triangles";
            }
        }
        return tris.empty() ? "" : "degenerate input produced triangles";
    }
    if (tris.empty())
    {
        return "non-degenerate input produced no triangles";
    }

    std::vector<DirEdge> es;
    es.reserve(3 * tris.size());
    std::vector<char> used(n, 0);
    Wide              area2 = Wide(0);
    for (size_t t = 0; t < tris.size(); ++t)
    {
        const Tri& tr = tris[t];
        for (int k = 0; k < 3; ++k)
        {
            if (tr[k] >= n) return "triangle index out of range";
            if (rep[tr[k]] != tr[k]) return "triangle references a duplicate (non-lowest) index";
            used[tr[k]] = 1;
        }
        if (tr[0] == tr[1] || tr[1] == tr[2] || tr[0] == tr[2]) return "triangle with repeated vertex";
        const Wide o = orient(pts[tr[0]], pts[tr[1]], pts[tr[2]]);
        if (sgn(o) <= 0) return "triangle not strictly counter-clockwise";
        area2 += o;
        for (int k = 0; k < 3; ++k)
        {
            const DirEdge e = { tr[k], tr[(k + 1) % 3], tr[(k + 2) % 3] };
            es.push_back(e);
        }
    }
    std::sort(es.begin(), es.end(), de_less);
    for (size_t i = 1; i < es.size(); ++i)
    {
        if (es[i].a == es[i - 1].a && es[i].b == es[i - 1].b)
        {
            return "directed edge used twice (overlapping triangles)";
        }
    }
    for (size_t i = 0; i < uniq.size(); ++i)
    {
        if (!used[uniq[i]]) return "distinct point not referenced by any triangle";
    }
    for (size_t i = 0; i < n; ++i)
    {
        if (qdt_point_is_used(q, i) != (used[i] ? 1u : 0u)) return "qdt_point_is_used disagrees with the triangles";
    }

    const std::vector<P> h = hull(up);
    Wide                 harea2 = Wide(0);
    for (size_t i = 0; i < h.size(); ++i)
    {
        const P& a = h[i];
        const P& b = h[(i + 1) % h.size()];
        harea2 += Wide(a.x) * Wide(b.y) - Wide(b.x) * Wide(a.y);
    }
    if (sgn(harea2) < 0)
    {
        harea2 = -harea2;
    }
    if (sgn(area2 - harea2) != 0)
    {
        return "triangle areas do not add up to the area of the convex hull";
    }

    for (size_t i = 0; i < es.size(); ++i)
    {
        const DirEdge& e = es[i];
        const DirEdge* rv = find_de(es, e.b, e.a);
        if (!rv)
        {
            /* boundary edge: every hull vertex must lie on its left (or on it) */
            for (size_t k = 0; k < h.size(); ++k)
            {
                if (sgn(orient(pts[e.a], pts[e.b], h[k])) < 0) return "boundary edge is not on the convex hull";
            }
        }
        else if (e.a < e.b)
        {
            if (!cons.empty() && cons.count(Edge(e.a, e.b))) continue;
            if (sgn(incircle(pts[e.a], pts[e.b], pts[e.c], pts[rv->c])) > 0) return "interior edge is not locally Delaunay";
        }
    }
    for (ConSet::const_iterator it = cons.begin(); it != cons.end(); ++it)
    {
        if (!find_de(es, it->first, it->second) && !find_de(es, it->second, it->first))
        {
            return "constrained edge missing from the triangulation";
        }
    }

    if (global && cons.empty())
    {
        for (size_t t = 0; t < tris.size(); ++t)
        {
            const Tri& tr = tris[t];
            for (size_t k = 0; k < uniq.size(); ++k)
            {
                const size_t i = uniq[k];
                if (i != tr[0] && i != tr[1] && i != tr[2] && sgn(incircle(pts[tr[0]], pts[tr[1]], pts[tr[2]], pts[i])) > 0)
                {
                    return "empty-circle property violated";
                }
            }
        }
    }
    return "";
}

/* ------------------------------------------------------------------------ */
/* Scenarios, part 1: API contract, inputs, exactness                       */
/* ------------------------------------------------------------------------ */

static TriList tris_of(const qdt_t* q)
{
    std::string e;
    return get_tris(q, &e);
}

static TriList mk(const size_t(*t)[3], size_t n)
{
    TriList r(n);
    for (size_t i = 0; i < n; ++i)
    {
        r[i][0] = t[i][0];
        r[i][1] = t[i][1];
        r[i][2] = t[i][2];
    }
    return r;
}

static void sc_examples(Scenario& r, const Config&)
{
    /* The quadrilateral of "Constrained triangulation". */
    {
        static const size_t before[2][3] = { { 3, 0, 1 },{ 1, 2, 3 } };
        static const size_t after[2][3] = { { 2, 0, 1 },{ 0, 2, 3 } };
        const qdt_vec2d     pts[4] = { { 0.0, 0.0 },{ 4.0, 0.0 },{ 3.0, 1.0 },{ 1.0, 1.0 } };
        qdt_t* q = qdt_create();
        r.expect(q != NULL, "qdt_create returned NULL");
        qdt_execute(q, pts, 4);
        r.expect(tri_sorted(tris_of(q)) == tri_sorted(mk(before, 2)), "quad: triangles before the constraint differ from the README");
        r.expect(qdt_insert_constraint(q, 0, 2) == 1, "quad: constraint 0-2 was not accepted");
        r.expect(tri_sorted(tris_of(q)) == tri_sorted(mk(after, 2)), "quad: triangles after the constraint differ from the README");
        qdt_destroy(q);
    }

    /* The six-point program of "Example". */
    {
        static const size_t listing[6][3] = { { 1, 3, 4 },{ 0, 1, 4 },{ 2, 0, 5 },{ 0, 4, 5 },{ 4, 3, 5 },{ 3, 2, 5 } };
        const qdt_vec2d     pts[6] = { { 0.0, 0.0 },{ 1.0, 0.0 },{ 0.0, 1.0 },{ 1.0, 1.0 },{ 0.5, 0.5 },{ 0.25, 0.75 } };
        for (int par = 0; par < 2; ++par)
        {
            qdt_t* q = qdt_create();
            exec_fn(par != 0)(q, pts, 4);
            ins_fn(par != 0)(q, pts[4]);
            ins_fn(par != 0)(q, pts[5]);
            r.expect(qdt_triangles_size(q) == 6, "six-point example: expected 6 triangles");
            const TriList got = tris_of(q);
            r.expect(tri_sorted(got) == tri_sorted(mk(listing, 6)), "six-point example: triangles differ from the README listing");
            if (par == 0)
            {
                bool same = got.size() == 6;
                for (size_t i = 0; same && i < 6; ++i)
                {
                    same = got[i][0] == listing[i][0] && got[i][1] == listing[i][1] && got[i][2] == listing[i][2];
                }
                r.note = std::string("The listing of the six-point example is ") +
                    (same ? "reproduced exactly, in the same order." : "reproduced as a set; the order differs from the README (documented as unspecified).");
            }
            qdt_destroy(q);
        }
    }

    /* The input array is copied, not kept and not modified. */
    for (int par = 0; par < 2; ++par)
    {
        std::vector<qdt_vec2d> d;
        d.push_back(qdt_vec2d{ 0.0, 0.0 });
        d.push_back(qdt_vec2d{ 4.0, 0.0 });
        d.push_back(qdt_vec2d{ 3.0, 1.0 });
        d.push_back(qdt_vec2d{ 1.0, 1.0 });
        const std::vector<qdt_vec2d> keep = d;
        qdt_t* q = qdt_create();
        exec_fn(par != 0)(q, d.data(), d.size());
        bool unchanged = true;
        for (size_t i = 0; i < d.size(); ++i)
        {
            unchanged = unchanged && d[i].x == keep[i].x && d[i].y == keep[i].y;
        }
        r.expect(unchanged, "a batch call modified the input array");
        for (size_t i = 0; i < d.size(); ++i)
        {
            d[i].x = -999.0;
            d[i].y = -999.0;
        }
        d.clear();
        d.shrink_to_fit();
        bool ok = qdt_triangles_size(q) == 2;
        for (size_t i = 0; i < keep.size(); ++i)
        {
            qdt_vec2d g = { 0.0, 0.0 };
            ok = ok && qdt_get_point(q, &g, i) == 1 && g.x == keep[i].x && g.y == keep[i].y;
        }
        r.expect(ok, "the context depends on the input array after the call");
        qdt_destroy(q);
    }
}

static void sc_api(Scenario& r, const Config&)
{
    /* NULL context */
    r.expect(qdt_insert_point(NULL, qdt_vec2d{ 0, 0 }) == 0, "insert_point(NULL) != 0");
    r.expect(qdt_insert_point_parallel(NULL, qdt_vec2d{ 0, 0 }) == 0, "insert_point_parallel(NULL) != 0");
    r.expect(qdt_insert_constraint(NULL, 0, 1) == 0, "insert_constraint(NULL) != 0");
    const qdt_edge ed = { 0, 1 };
    r.expect(qdt_insert_constraints_batch(NULL, &ed, 1) == 0, "insert_constraints_batch(NULL) != 0");
    r.expect(qdt_triangles_size(NULL) == 0, "triangles_size(NULL) != 0");
    size_t t[3] = { 7, 8, 9 };
    r.expect(qdt_get_triangle(NULL, t, 0) == 0 && t[0] == 7, "get_triangle(NULL) misbehaves");
    qdt_execute(NULL, NULL, 0);
    qdt_execute_parallel(NULL, NULL, 0);
    qdt_destroy(NULL);
    {
        qdt_vec2d tmp = { 5.0, 6.0 };
        r.expect(qdt_points_size(NULL) == 0, "points_size(NULL) != 0");
        r.expect(qdt_get_point(NULL, &tmp, 0) == 0 && tmp.x == 5.0 && tmp.y == 6.0, "get_point(NULL) misbehaves");
        r.expect(qdt_point_is_used(NULL, 0) == 0, "point_is_used(NULL) != 0");
    }

    /* fresh context */
    qdt_t* q = qdt_create();
    r.expect(qdt_triangles_size(q) == 0 && qdt_points_size(q) == 0, "a fresh context is not empty");
    r.expect(qdt_insert_constraint(q, 0, 1) == 0, "constraint accepted on an empty context");
    r.expect(qdt_insert_constraints_batch(q, &ed, 1) == 0, "batch constraint accepted on an empty context");
    r.expect(qdt_insert_constraints_batch(q, NULL, 0) == 0, "insert_constraints_batch(NULL, 0) != 0");

    /* queries on a const context */
    const qdt_vec2d sq[4] = { { 0, 0 },{ 4, 0 },{ 3, 1 },{ 1, 1 } };
    qdt_execute(q, sq, 4);
    {
        const qdt_t* cq = q;
        qdt_vec2d    got = { -1.0, -1.0 };
        size_t       tt[3] = { 7, 8, 9 };
        r.expect(qdt_triangles_size(cq) == 2, "triangles_size(const) wrong");
        r.expect(qdt_get_triangle(cq, tt, 0) == 1, "get_triangle(const) failed");
        r.expect(qdt_get_triangle(cq, NULL, 0) == 0, "get_triangle(NULL output) != 0");
        r.expect(qdt_points_size(cq) == 4, "points_size after execute != 4");
        r.expect(qdt_get_point(cq, &got, 2) == 1 && got.x == 3.0 && got.y == 1.0, "get_point(2) wrong");
        got.x = -1.0;
        got.y = -1.0;
        r.expect(qdt_get_point(cq, &got, 4) == 0 && got.x == -1.0 && got.y == -1.0, "get_point(out of range) misbehaves");
        r.expect(qdt_get_point(cq, NULL, 0) == 0, "get_point(NULL output) != 0");
        bool all = true;
        for (size_t i = 0; i < 4; ++i)
        {
            all = all && qdt_point_is_used(cq, i) == 1;
        }
        r.expect(all, "a vertex of the quadrilateral is reported as unused");
        r.expect(qdt_point_is_used(cq, 4) == 0, "point_is_used(out of range) != 0");
    }

    /* duplicates and degenerate contexts in the point queries */
    {
        const qdt_vec2d dup[5] = { { 0, 0 },{ 1, 0 },{ 0, 0 },{ 0, 1 },{ 1, 0 } };
        const size_t    want[5] = { 1, 1, 0, 1, 0 };
        qdt_t* g = qdt_create();
        qdt_execute(g, dup, 5);
        r.expect(qdt_points_size(g) == 5, "duplicates: points_size != 5");
        bool ok = true;
        for (size_t i = 0; i < 5; ++i)
        {
            ok = ok && qdt_point_is_used(g, i) == want[i];
        }
        r.expect(ok, "duplicates: point_is_used wrong");
        r.expect(qdt_insert_point(g, qdt_vec2d{ 0.0, 1.0 }) == 6, "duplicates: inserting a duplicate returned the wrong count");
        r.expect(qdt_point_is_used(g, 5) == 0 && qdt_point_is_used(g, 3) == 1, "duplicates: inserted duplicate reported as used");
        r.expect(qdt_insert_point(g, qdt_vec2d{ 2.0, 2.0 }) == 7, "duplicates: insert returned the wrong count");
        r.expect(qdt_point_is_used(g, 6) == 1, "duplicates: inserted point reported as unused");
        qdt_execute(g, dup, 2);
        r.expect(qdt_points_size(g) == 2 && qdt_point_is_used(g, 0) == 0 && qdt_point_is_used(g, 1) == 0,
            "degenerate: point queries wrong");
        qdt_destroy(g);
    }

    /* constraint end points: points inserted after the batch call, their duplicates, later points */
    {
        const qdt_vec2d b[4] = { { 0, 0 },{ 4, 0 },{ 4, 4 },{ 0, 4 } };
        qdt_t* g = qdt_create();
        qdt_execute(g, b, 4);
        const size_t i4 = qdt_insert_point(g, qdt_vec2d{ 2.0, 1.0 }) - 1;
        const size_t i5 = qdt_insert_point(g, qdt_vec2d{ 2.0, 1.0 }) - 1; /* duplicate of 4 */
        const size_t i6 = qdt_insert_point(g, qdt_vec2d{ 1.0, 3.0 }) - 1;
        r.expect(i4 == 4 && i5 == 5 && i6 == 6, "constraint lookup: unexpected point indices");
        r.expect(qdt_insert_constraint(g, i4, i6) == 1, "constraint between two inserted points rejected");
        r.expect(qdt_insert_constraint(g, i5, 0) == 1, "constraint through a duplicate's index rejected");
        r.expect(qdt_insert_constraint(g, i5, i4) == 0, "constraint between duplicates accepted");
        const qdt_edge be[2] = { { i5, i6 },{ i4, 2 } };
        r.expect(qdt_insert_constraints_batch(g, be, 2) == 2, "batch with inserted points: wrong count");
        const size_t i7 = qdt_insert_point(g, qdt_vec2d{ 0.0, 2.0 }) - 1;
        r.expect(i7 == 7, "constraint lookup: unexpected index of the later point");
        r.expect(qdt_insert_constraint(g, i7, i6) == 1, "constraint to a point inserted after earlier constraints rejected");
        qdt_destroy(g);
    }

    /* constraint rules and replacement by the next batch call */
    r.expect(qdt_insert_constraint(q, 0, 99) == 0, "constraint with an invalid index accepted");
    r.expect(qdt_insert_constraint(q, 99, 98) == 0, "constraint with two invalid indices accepted");
    r.expect(qdt_insert_constraint(q, static_cast<size_t>(-1), 0) == 0, "constraint with index SIZE_MAX accepted");
    r.expect(qdt_insert_constraint(q, 2, 2) == 0, "constraint with v0 == v1 accepted");
    r.expect(qdt_insert_constraint(q, 0, 2) == 1, "constraint 0-2 rejected");
    r.expect(qdt_insert_constraint(q, 2, 0) == 1, "inserting an existing constraint (reversed) rejected");
    r.expect(qdt_insert_constraint(q, 1, 3) == 0, "constraint crossing an existing constraint accepted");
    qdt_execute(q, sq, 4);
    r.expect(qdt_insert_constraint(q, 1, 3) == 1, "constraints were not discarded by qdt_execute");
    qdt_execute(q, NULL, 0);
    r.expect(qdt_triangles_size(q) == 0 && qdt_points_size(q) == 0, "execute with n = 0 left data behind");
    for (int i = 0; i < 5; ++i)
    {
        qdt_insert_point(q, qdt_vec2d{ static_cast<double>(i), static_cast<double>(i * i) });
    }
    r.expect(qdt_triangles_size(q) == 3, "insertions into an emptied context: wrong triangle count");
    qdt_destroy(q);
}

static void sc_guards(Scenario& r, const Config&)
{
    qdt_vec2d ok4[4] = { { 0, 0 },{ 4, 0 },{ 3, 1 },{ 1, 1 } };
    qdt_vec2d nan4[4] = { { 0, 0 },{ 4, 0 },{ std::nan(""), 1 },{ 1, 1 } };
    qdt_vec2d inf4[4] = { { 0, 0 },{ 4, 0 },{ 3, 1 },{ HUGE_VAL, 1 } };
    qdt_vec2d ninf4[4] = { { 0, 0 },{ 4, 0 },{ 3, 1 },{ 1, -HUGE_VAL } };
    for (int par = 0; par < 2; ++par)
    {
        const ExecFn run = exec_fn(par != 0);
        qdt_t* g = qdt_create();
        run(g, ok4, 4);
        r.expect(qdt_triangles_size(g) == 2, "guards: setup did not give 2 triangles");
        run(g, NULL, 4);
        r.expect(qdt_triangles_size(g) == 0 && qdt_points_size(g) == 0, "a NULL array with n > 0 did not empty the context");
        run(g, ok4, 4);
        run(g, nan4, 4);
        r.expect(qdt_triangles_size(g) == 0 && qdt_points_size(g) == 0, "NaN input did not empty the context");
        run(g, ok4, 4);
        run(g, inf4, 4);
        r.expect(qdt_triangles_size(g) == 0 && qdt_points_size(g) == 0, "+infinity input did not empty the context");
        run(g, ok4, 4);
        run(g, ninf4, 4);
        r.expect(qdt_triangles_size(g) == 0 && qdt_points_size(g) == 0, "-infinity input did not empty the context");
        run(g, ok4, 4);
        run(g, ok4, static_cast<size_t>(100000001)); /* rejected by the size limit before any point is read */
        r.expect(qdt_triangles_size(g) == 0 && qdt_points_size(g) == 0, "n > 100,000,000 did not empty the context");
        run(g, ok4, 4);
        r.expect(qdt_triangles_size(g) == 2, "the context is unusable after a rejected call");

        const InsFn ins = ins_fn(par != 0);
        r.expect(ins(g, qdt_vec2d{ std::nan(""), 0.0 }) == 4, "a NaN point was inserted");
        r.expect(ins(g, qdt_vec2d{ 0.0, HUGE_VAL }) == 4, "an infinite point was inserted");
        r.expect(qdt_triangles_size(g) == 2 && qdt_points_size(g) == 4, "a rejected point changed the context");
        r.expect(ins(g, qdt_vec2d{ 2.0, 0.5 }) == 5, "a valid point after rejected ones: wrong count");
        qdt_destroy(g);
    }
}

static void sc_structured(Scenario& r, const Config&)
{
    struct Case
    {
        const char* name;
        std::vector<P> pts;
    };
    std::vector<Case> cases;
    const double      kPi = 3.14159265358979323846;
    {
        std::vector<P> p;
        for (int i = 0; i < 40; ++i)
            for (int j = 0; j < 40; ++j) p.push_back(P{ i, j });
        cases.push_back(Case{ "grid 40x40", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 50; ++i) p.push_back(P{ i, 2 * i });
        cases.push_back(Case{ "collinear diagonal", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 50; ++i) p.push_back(P{ i, 7 });
        cases.push_back(Case{ "collinear horizontal", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 50; ++i) p.push_back(P{ 7, i });
        cases.push_back(Case{ "collinear vertical", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 50; ++i) p.push_back(P{ i, 7 });
        p.push_back(P{ 25, 8 });
        cases.push_back(Case{ "line plus one point", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 30; ++i) p.push_back(P{ i, static_cast<ll>(i) * i });
        cases.push_back(Case{ "parabola", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 100; ++i) p.push_back(P{ 5, 5 });
        cases.push_back(Case{ "all duplicates", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 200; ++i)
        {
            const double a = 2 * kPi * i / 200;
            p.push_back(P{ static_cast<ll>(std::llround((1 << 27) * std::cos(a))), static_cast<ll>(std::llround((1 << 27) * std::sin(a))) });
        }
        cases.push_back(Case{ "near circle, radius 2^27", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 100; ++i)
        {
            p.push_back(P{ i, 0 });
            p.push_back(P{ i, 1 });
        }
        cases.push_back(Case{ "two rows", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 20; ++i)
        {
            p.push_back(P{ 0, 0 });
            p.push_back(P{ 1, 0 });
            p.push_back(P{ 0, 1 });
        }
        cases.push_back(Case{ "triangle with duplicates", p });
    }
    cases.push_back(Case{ "empty", std::vector<P>() });
    cases.push_back(Case{ "one point", std::vector<P>(1, P{ 1, 1 }) });
    {
        std::vector<P> p;
        p.push_back(P{ 1, 1 });
        p.push_back(P{ 2, 3 });
        cases.push_back(Case{ "two points", p });
    }
    {
        std::vector<P> p;
        for (int i = 0; i < 3; ++i) p.push_back(P{ i, i });
        cases.push_back(Case{ "three collinear points", p });
    }
    {
        std::vector<P> p;
        p.push_back(P{ 0, 0 });
        p.push_back(P{ 5, 0 });
        p.push_back(P{ 0, 5 });
        cases.push_back(Case{ "three points", p });
    }
    {
        const ll       m = 1 << 28;
        std::vector<P> p;
        p.push_back(P{ -m, -m });
        p.push_back(P{ m, -m });
        p.push_back(P{ 0, m });
        p.push_back(P{ 1, 1 });
        p.push_back(P{ -1, 2 });
        cases.push_back(Case{ "coordinates up to 2^28", p });
    }

    for (int par = 0; par < 2; ++par)
    {
        for (size_t ci = 0; ci < cases.size(); ++ci)
        {
            for (int shuffle = 0; shuffle < 5; ++shuffle)
            {
                std::vector<P> pts = cases[ci].pts;
                std::mt19937_64 rng(77 + shuffle);
                if (shuffle)
                {
                    std::shuffle(pts.begin(), pts.end(), rng);
                }
                const std::vector<qdt_vec2d> d = to_d(pts);
                qdt_t* q = qdt_create();
                exec_fn(par != 0)(q, d.empty() ? NULL : d.data(), d.size());
                std::string e = verify(q, pts, ConSet(), pts.size() < 1000);
                if (!e.empty())
                {
                    e = std::string(cases[ci].name) + ": " + e + (par ? " (parallel)" : "");
                }
                r.finish(static_cast<uint64_t>(shuffle), e);
                qdt_destroy(q);
            }
        }
    }
}

static void sc_batch(Scenario& r, const Config& c)
{
    static const ll ranges[] = { 2, 3, 4, 6, 10, 30, 100, 10000, 1 << 20 };
    for (int par = 0; par < 2; ++par)
    {
        for (uint64_t s = 0; s < c.iters; ++s)
        {
            const uint64_t  sd = seed_of(c, 1000 + s);
            std::mt19937_64 rng(sd);
            const size_t    n = rng() % 70;
            const ll        R = ranges[rng() % (sizeof(ranges) / sizeof(ranges[0]))];
            const std::vector<P>         pts = gen(rng, n, R);
            const std::vector<qdt_vec2d> d = to_d(pts);
            qdt_t* q = qdt_create();
            exec_fn(par != 0)(q, d.data(), n);
            std::string e = verify(q, pts, ConSet(), true);
            if (!e.empty())
            {
                e += fmt(" (n=%zu, range=%lld%s)", n, R, par ? ", parallel" : "");
            }
            r.finish(sd, e);
            qdt_destroy(q);
        }
    }
}

/* One set that mixes coordinates of very different magnitude. */
static void sc_magnitudes(Scenario& r, const Config& c)
{
    static const int levels[] = { 0, 2, 5, 9, 14, 20, 26, 28 };
    for (uint64_t s = 0; s < c.iters; ++s)
    {
        const uint64_t  sd = seed_of(c, 41000 + s);
        std::mt19937_64 rng(sd);
        const int       mode = static_cast<int>(s % 3);
        const size_t    n = 4 + rng() % 60;
        const ll        offx = ((rng() & 1) ? 1 : -1) * (static_cast<ll>(1) << 27);
        const ll        offy = ((rng() & 1) ? 1 : -1) * (static_cast<ll>(1) << 27);
        std::vector<P>  pts;
        for (size_t i = 0; i < n; ++i)
        {
            const bool far_cluster = mode == 1 || (mode == 2 && (i & 1));
            if (far_cluster)
            {
                const ll b = static_cast<ll>(1) << (rng() % 9);
                pts.push_back(P{ offx + signed_rand(rng, b), offy + signed_rand(rng, b) });
            }
            else
            {
                const ll b = static_cast<ll>(1) << levels[rng() % (sizeof(levels) / sizeof(levels[0]))];
                pts.push_back(P{ signed_rand(rng, b), signed_rand(rng, b) });
            }
        }
        const std::vector<qdt_vec2d> d = to_d(pts);
        qdt_t* q = qdt_create();
        exec_fn((s & 1) != 0)(q, d.data(), n);
        std::string e = verify(q, pts, ConSet(), true);
        if (!e.empty())
        {
            e += fmt(" (mode=%d, n=%zu)", mode, n);
        }
        r.finish(sd, e);
        qdt_destroy(q);
    }
}

/* Exactness under overall scaling: coordinates times 2^k are exact doubles, so
   the integer verifier applies to the original integers. */
static void sc_scale(Scenario& r, const Config& c)
{
    static const int ks[] = { -1000, -600, -300, -64, 64, 300, 600, 900 };
    const uint64_t   per = c.quick ? 20 : 150;
    for (size_t ki = 0; ki < sizeof(ks) / sizeof(ks[0]); ++ki)
    {
        for (uint64_t s = 0; s < per; ++s)
        {
            const uint64_t  sd = seed_of(c, 31000 + s);
            std::mt19937_64 rng(sd);
            const ll        R = (rng() % 2) ? 6 : 1000;
            const size_t    n = 3 + rng() % 40;
            const std::vector<P>         pts = gen(rng, n, R);
            const std::vector<qdt_vec2d> d = to_d(pts, std::ldexp(1.0, ks[ki]));
            qdt_t* q = qdt_create();
            exec_fn((s & 1) != 0)(q, d.data(), n);
            std::string e = verify(q, pts, ConSet(), true, std::ldexp(1.0, ks[ki]));
            if (!e.empty())
            {
                e += fmt(" (scale 2^%d)", ks[ki]);
            }
            r.finish(sd, e);
            qdt_destroy(q);
        }
    }
}

static void sc_circle(Scenario& r, const Config& c)
{
    const double   kPi = 3.14159265358979323846;
    const uint64_t cnt = c.quick ? 100 : 400;
    for (uint64_t s = 0; s < cnt; ++s)
    {
        const uint64_t  sd = seed_of(c, 55000 + s);
        std::mt19937_64 rng(sd);
        const size_t    n = 8 + rng() % 120;
        const ll        Rr = static_cast<ll>(1) << (20 + rng() % 8);
        std::vector<P>  pts;
        for (size_t i = 0; i < n; ++i)
        {
            const double a = 2 * kPi * static_cast<double>(rng() % 100000) / 100000.0;
            pts.push_back(P{ static_cast<ll>(std::llround(Rr * std::cos(a))), static_cast<ll>(std::llround(Rr * std::sin(a))) });
        }
        const std::vector<qdt_vec2d> d = to_d(pts);
        qdt_t* q = qdt_create();
        qdt_execute(q, d.data(), n);
        std::string e = verify(q, pts, ConSet(), true);
        for (int i = 0; i < 3 && e.empty(); ++i) /* a few insertions on top */
        {
            const double a = 2 * kPi * static_cast<double>(rng() % 100000) / 100000.0;
            const P      p = { static_cast<ll>(std::llround(Rr * std::cos(a))), static_cast<ll>(std::llround(Rr * std::sin(a))) };
            qdt_insert_point(q, qdt_vec2d{ static_cast<double>(p.x), static_cast<double>(p.y) });
            pts.push_back(p);
            e = verify(q, pts, ConSet(), true);
            if (!e.empty())
            {
                e += " (after an insertion)";
            }
        }
        r.finish(sd, e);
        qdt_destroy(q);
    }
}

/* ------------------------------------------------------------------------ */
/* Scenarios, part 2: insertion and constraints                             */
/* ------------------------------------------------------------------------ */

static void sc_insert(Scenario& r, const Config& c)
{
    static const ll ranges[] = { 2, 3, 5, 8, 20, 100, 10000 };
    for (int par = 0; par < 2; ++par)
    {
        const InsFn ins = ins_fn(par != 0);
        for (uint64_t s = 0; s < c.iters; ++s)
        {
            const uint64_t  sd = seed_of(c, 5000 + s);
            std::mt19937_64 rng(sd);
            const ll        R = ranges[rng() % (sizeof(ranges) / sizeof(ranges[0]))];
            const size_t    n0 = rng() % 40;
            const size_t    k = 1 + rng() % 25;
            std::vector<P>  pts = gen(rng, n0, R);
            qdt_t* q = qdt_create();
            const bool      do_exec = (rng() % 8) != 0; /* sometimes insert into a context that was never filled */
            if (do_exec)
            {
                const std::vector<qdt_vec2d> d = to_d(pts);
                exec_fn(par != 0)(q, d.data(), n0);
            }
            else
            {
                pts.clear();
            }
            std::string msg;
            for (size_t i = 0; i < k && msg.empty(); ++i)
            {
                P p;
                if (!pts.empty() && rng() % 5 == 0)
                {
                    p = pts[rng() % pts.size()]; /* duplicate */
                }
                else if (pts.size() >= 2 && rng() % 4 == 0)
                {
                    /* lattice point on the segment between two points, when there is one */
                    const P a = pts[rng() % pts.size()], b = pts[rng() % pts.size()];
                    p = P{ (a.x + b.x) / 2, (a.y + b.y) / 2 };
                    if ((a.x + b.x) % 2 || (a.y + b.y) % 2) p = gen(rng, 1, R)[0];
                }
                else
                {
                    p = gen(rng, 1, R)[0];
                }
                const size_t ret = ins(q, qdt_vec2d{ static_cast<double>(p.x), static_cast<double>(p.y) });
                pts.push_back(p);
                if (ret != pts.size())
                {
                    msg = fmt("insert_point returned %zu, expected %zu", ret, pts.size());
                }
                else
                {
                    msg = verify(q, pts, ConSet(), true);
                }
                if (!msg.empty())
                {
                    msg += fmt(" (after insertion #%zu, n0=%zu, range=%lld, %s, %s)", i, n0, R, do_exec ? "batch first" : "empty context",
                        par ? "parallel" : "serial");
                }
            }
            r.finish(sd, msg);
            qdt_destroy(q);
        }
    }
}

/* Insertions into a larger triangulation (the point location has to travel):
   duplicates of existing points, points on edges, points outside the hull. */
static void sc_insertlarge(Scenario& r, const Config& c)
{
    const uint64_t cnt = c.iters / 100 + 1;
    for (uint64_t s = 0; s < cnt; ++s)
    {
        const uint64_t  sd = seed_of(c, 77000 + s);
        std::mt19937_64 rng(sd);
        const ll        R = static_cast<ll>(1) << 20;
        const size_t    n0 = 2000 + rng() % 4000;
        std::vector<P>  pts = gen(rng, n0, R);
        const std::vector<qdt_vec2d> d = to_d(pts);
        qdt_t* q = qdt_create();
        qdt_execute(q, d.data(), n0);
        std::string msg;
        for (size_t i = 0; i < 400 && msg.empty(); ++i)
        {
            P p;
            switch (rng() % 4)
            {
            case 0: p = pts[rng() % pts.size()]; break;
            case 1:
            {
                const P a = pts[rng() % pts.size()], b = pts[rng() % pts.size()];
                p = P{ (a.x + b.x) / 2, (a.y + b.y) / 2 };
                if ((a.x + b.x) % 2 || (a.y + b.y) % 2) p = gen(rng, 1, R)[0];
                break;
            }
            case 2: p = P{ signed_rand(rng, R) , signed_rand(rng, R) }; break; /* may be outside the hull */
            default: p = gen(rng, 1, R)[0];
            }
            const size_t ret = qdt_insert_point(q, qdt_vec2d{ static_cast<double>(p.x), static_cast<double>(p.y) });
            pts.push_back(p);
            if (ret != pts.size())
            {
                msg = "insert_point returned the wrong count";
                break;
            }
            if (i % 100 == 99)
            {
                msg = verify(q, pts, ConSet(), false);
            }
        }
        if (msg.empty())
        {
            msg = verify(q, pts, ConSet(), false);
        }
        r.finish(sd, msg);
        qdt_destroy(q);
    }
}

/* ---- constraints ---- */

static std::vector<size_t> lowest_dup(const std::vector<P>& pts)
{
    const size_t        n = pts.size();
    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), static_cast<size_t>(0));
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (!coord_eq(pts[a], pts[b])) return coord_less(pts[a], pts[b]);
        return a < b;
        });
    std::vector<size_t> rep(n);
    for (size_t i = 0; i < n; ++i)
    {
        rep[order[i]] = (i > 0 && coord_eq(pts[order[i]], pts[order[i - 1]])) ? rep[order[i - 1]] : order[i];
    }
    return rep;
}

static Edge min_max(size_t a, size_t b)
{
    return Edge(std::min(a, b), std::max(a, b));
}

/* Whether the documented rules say that a constraint between the points a and b
   is accepted, given the constraints accepted so far (as lowest-duplicate pairs). */
static bool constraint_ok(const std::vector<P>& pts, const std::vector<size_t>& rep, const std::vector<Edge>& acc, size_t a,
    size_t b, bool have_tris)
{
    const size_t ra = rep[a], rb = rep[b];
    if (ra == rb || !have_tris)
    {
        return false;
    }
    for (size_t j = 0; j < pts.size(); ++j)
    {
        if (rep[j] == j && on_open_segment(pts[ra], pts[rb], pts[j]))
        {
            return false; /* passes exactly through a third point */
        }
    }
    for (size_t i = 0; i < acc.size(); ++i)
    {
        if (proper_cross(pts[ra], pts[rb], pts[acc[i].first], pts[acc[i].second]))
        {
            return false; /* crosses an accepted constraint */
        }
    }
    return true;
}

static ll gcd_ll(ll a, ll b)
{
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;
    while (b)
    {
        const ll t = a % b;
        a = b;
        b = t;
    }
    return a;
}

/* Random constraints against an independent oracle for the documented rules;
   with do_insert, point insertions are interleaved (constraints are kept and
   split by points inserted on them). Runs the single-call and the batch API. */
static void sc_constraints(Scenario& r, const Config& c, bool do_insert)
{
    static const ll ranges[] = { 4, 6, 10, 20, 50, 1000 };
    for (int batch = 0; batch < 2; ++batch)
    {
        for (uint64_t s = 0; s < c.iters; ++s)
        {
            const uint64_t  sd = seed_of(c, 9000 + s + (batch ? 100000 : 0) + (do_insert ? 200000 : 0));
            std::mt19937_64 rng(sd);
            const ll        R = ranges[rng() % (sizeof(ranges) / sizeof(ranges[0]))];
            const size_t    n = 4 + rng() % 40;
            std::vector<P>  pts = gen(rng, n, R);
            const std::vector<qdt_vec2d> d = to_d(pts);
            qdt_t* q = qdt_create();
            qdt_execute(q, d.data(), n);
            std::vector<Edge> acc; /* accepted constraints as lowest-duplicate pairs */
            std::string       msg;
            const int         rounds = 1 + static_cast<int>(rng() % 6);
            for (int rd = 0; rd < rounds && msg.empty(); ++rd)
            {
                const size_t              m = 1 + rng() % 8;
                const std::vector<size_t> rep = lowest_dup(pts);
                const bool                tris = has_tris(pts);
                std::vector<qdt_edge>     edges;
                std::vector<char>         expect;
                std::vector<Edge>         cand;
                std::vector<Edge>         acc_local = acc;
                for (size_t i = 0; i < m; ++i)
                {
                    const size_t a = rng() % pts.size(), b = rng() % pts.size();
                    const bool   ok = constraint_ok(pts, rep, acc_local, a, b, tris);
                    if (ok)
                    {
                        acc_local.push_back(min_max(rep[a], rep[b]));
                    }
                    const qdt_edge e = { a, b };
                    edges.push_back(e);
                    expect.push_back(ok ? 1 : 0);
                    cand.push_back(Edge(rep[a], rep[b]));
                }
                if (batch)
                {
                    const size_t accepted = qdt_insert_constraints_batch(q, edges.data(), edges.size());
                    size_t       want = 0;
                    for (size_t i = 0; i < m; ++i) want += expect[i] ? 1 : 0;
                    if (accepted != want)
                    {
                        msg = fmt("batch accepted %zu constraints, expected %zu", accepted, want);
                    }
                }
                else
                {
                    for (size_t i = 0; i < m && msg.empty(); ++i)
                    {
                        const TriList before = tri_sorted(tris_of(q));
                        const size_t  rv = qdt_insert_constraint(q, edges[i].v0, edges[i].v1);
                        if (rv > 1)
                        {
                            msg = "insert_constraint returned a value above 1";
                        }
                        else if ((rv == 1) != (expect[i] != 0))
                        {
                            msg = fmt("%s constraint (%zu, %zu)", rv ? "unexpected accept of" : "wrongly rejected", edges[i].v0, edges[i].v1);
                        }
                        else if (rv == 0 && tri_sorted(tris_of(q)) != before)
                        {
                            msg = "a rejected constraint changed the triangulation";
                        }
                    }
                }
                if (!msg.empty())
                {
                    msg += fmt(" (%s, n=%zu, range=%lld)", batch ? "batch" : "single calls", pts.size(), R);
                    break;
                }
                for (size_t i = 0; i < m; ++i)
                {
                    if (expect[i]) acc.push_back(min_max(cand[i].first, cand[i].second));
                }
                msg = verify(q, pts, ConSet(acc.begin(), acc.end()), false);
                if (!msg.empty())
                {
                    msg += fmt(" (after round %d of constraints, %s, n=%zu, range=%lld)", rd, batch ? "batch" : "single calls", pts.size(), R);
                    break;
                }

                if (!do_insert)
                {
                    continue;
                }
                const size_t k = 1 + rng() % 4;
                for (size_t i = 0; i < k && msg.empty(); ++i)
                {
                    P p;
                    if (!acc.empty() && rng() % 2)
                    {
                        /* on a constraint, if there is a lattice point there */
                        const Edge e = acc[rng() % acc.size()];
                        const P    a = pts[e.first], b = pts[e.second];
                        const ll   dx = b.x - a.x, dy = b.y - a.y;
                        const ll   g = gcd_ll(dx, dy);
                        if (g > 1)
                        {
                            const ll t = 1 + static_cast<ll>(rng() % static_cast<uint64_t>(g - 1));
                            p = P{ a.x + dx / g * t, a.y + dy / g * t };
                        }
                        else
                        {
                            p = gen(rng, 1, R)[0];
                        }
                    }
                    else
                    {
                        p = gen(rng, 1, R)[0];
                    }
                    qdt_insert_point(q, qdt_vec2d{ static_cast<double>(p.x), static_cast<double>(p.y) });
                    const size_t idx = pts.size();
                    pts.push_back(p);
                    bool dup = false;
                    for (size_t j = 0; j < idx; ++j)
                    {
                        if (coord_eq(pts[j], p)) dup = true;
                    }
                    if (!dup)
                    {
                        /* constraints that p lies on are split in two */
                        std::vector<Edge> next;
                        for (size_t j = 0; j < acc.size(); ++j)
                        {
                            if (on_open_segment(pts[acc[j].first], pts[acc[j].second], p))
                            {
                                next.push_back(min_max(acc[j].first, idx));
                                next.push_back(min_max(acc[j].second, idx));
                            }
                            else
                            {
                                next.push_back(acc[j]);
                            }
                        }
                        acc = next;
                    }
                    msg = verify(q, pts, ConSet(acc.begin(), acc.end()), false);
                    if (!msg.empty())
                    {
                        msg += fmt(" (after the insertion of point %zu%s, n=%zu, range=%lld)", idx, dup ? ", a duplicate" : "", pts.size(), R);
                    }
                }
            }
            r.finish(sd, msg);
            qdt_destroy(q);
        }
    }
}

static void sc_cdt(Scenario& r, const Config& c)
{
    sc_constraints(r, c, false);
}

static void sc_cdtins(Scenario& r, const Config& c)
{
    sc_constraints(r, c, true);
}

/* No 3 collinear and no 4 cocircular points: the constrained Delaunay
   triangulation of the points and a compatible set of constraints is unique. */
static bool general_position(const std::vector<P>& u)
{
    const size_t n = u.size();
    for (size_t a = 0; a < n; ++a)
        for (size_t b = a + 1; b < n; ++b)
            for (size_t c = b + 1; c < n; ++c)
                if (sgn(orient(u[a], u[b], u[c])) == 0) return false;
    for (size_t a = 0; a < n; ++a)
        for (size_t b = a + 1; b < n; ++b)
            for (size_t c = b + 1; c < n; ++c)
                for (size_t d = c + 1; d < n; ++d)
                    if (sgn(incircle(u[a], u[b], u[c], u[d])) == 0) return false;
    return true;
}

static qdt_t* make_ctx(const std::vector<P>& pts)
{
    const std::vector<qdt_vec2d> d = to_d(pts);
    qdt_t* q = qdt_create();
    qdt_execute(q, d.data(), d.size());
    return q;
}

/* README: "Same as calling qdt_insert_constraint for each entry in order" and
   "the order of the constraints does not matter" in general position. */
static void sc_cdtequiv(Scenario& r, const Config& c)
{
    const uint64_t cnt = c.iters / 3 + 1;
    for (uint64_t s = 0; s < cnt; ++s)
    {
        const uint64_t  sd = seed_of(c, 61000 + s);
        std::mt19937_64 rng(sd);
        const size_t    n = 8 + rng() % 9;
        const std::vector<P> pts = gen(rng, n, static_cast<ll>(1) << 26);
        if (!general_position(pts))
        {
            continue; /* not a case: the comparison below needs a unique result */
        }
        const size_t          m = 3 + rng() % 8;
        std::vector<qdt_edge> list;
        for (size_t i = 0; i < m; ++i)
        {
            const qdt_edge e = { rng() % n, rng() % n };
            list.push_back(e);
        }
        std::string msg;

        /* 1. one batch call == single calls in the same order (rejections included) */
        {
            qdt_t* qs = make_ctx(pts);
            qdt_t* qb = make_ctx(pts);
            size_t  acc_s = 0;
            ConSet  cons;
            for (size_t i = 0; i < m; ++i)
            {
                if (qdt_insert_constraint(qs, list[i].v0, list[i].v1) == 1)
                {
                    ++acc_s;
                    cons.insert(min_max(list[i].v0, list[i].v1));
                }
            }
            const size_t acc_b = qdt_insert_constraints_batch(qb, list.data(), list.size());
            if (acc_s != acc_b)
            {
                msg = fmt("batch accepted %zu, single calls %zu", acc_b, acc_s);
            }
            else if (tri_sorted(tris_of(qs)) != tri_sorted(tris_of(qb)))
            {
                msg = "batch and single calls gave different triangulations";
            }
            else
            {
                msg = verify(qb, pts, cons, false);
            }
            qdt_destroy(qs);
            qdt_destroy(qb);
        }

        /* 2. a compatible set gives the same triangulation in every order */
        if (msg.empty())
        {
            std::vector<size_t> rep(n);
            std::iota(rep.begin(), rep.end(), static_cast<size_t>(0));
            std::vector<Edge> K;
            for (size_t i = 0; i < m; ++i)
            {
                const Edge e = min_max(list[i].v0, list[i].v1);
                if (std::find(K.begin(), K.end(), e) == K.end() && constraint_ok(pts, rep, K, e.first, e.second, true))
                {
                    K.push_back(e);
                }
            }
            if (K.size() >= 2)
            {
                std::vector<qdt_edge> fwd, rev, shuf;
                for (size_t i = 0; i < K.size(); ++i)
                {
                    const qdt_edge e = { K[i].first, K[i].second };
                    fwd.push_back(e);
                }
                rev.assign(fwd.rbegin(), fwd.rend());
                shuf = fwd;
                std::shuffle(shuf.begin(), shuf.end(), rng);

                qdt_t* qf = make_ctx(pts);
                qdt_t* qr = make_ctx(pts);
                qdt_t* qh = make_ctx(pts);
                qdt_t* qg = make_ctx(pts);
                bool    all = true;
                for (size_t i = 0; i < fwd.size(); ++i)
                {
                    all = all && qdt_insert_constraint(qf, fwd[i].v0, fwd[i].v1) == 1;
                    all = all && qdt_insert_constraint(qr, rev[i].v0, rev[i].v1) == 1;
                }
                all = all && qdt_insert_constraints_batch(qh, shuf.data(), shuf.size()) == shuf.size();
                all = all && qdt_insert_constraints_batch(qg, fwd.data(), fwd.size()) == fwd.size();
                if (!all)
                {
                    msg = "a compatible constraint was rejected";
                }
                else
                {
                    const TriList ref = tri_sorted(tris_of(qf));
                    if (tri_sorted(tris_of(qr)) != ref || tri_sorted(tris_of(qh)) != ref || tri_sorted(tris_of(qg)) != ref)
                    {
                        msg = "the order of the constraints changed the triangulation";
                    }
                    else
                    {
                        msg = verify(qf, pts, ConSet(K.begin(), K.end()), false);
                    }
                }
                qdt_destroy(qf);
                qdt_destroy(qr);
                qdt_destroy(qh);
                qdt_destroy(qg);
            }
        }
        r.finish(sd, msg);
    }
}

/* Base set large enough that the first insertion builds the internal mesh in
   parallel (at least 50,000 triangles), then insertions and constraints that
   use the point -> vertex table built there. */
static void parmesh_case(Scenario& r, bool par_exec, bool par_ins, uint64_t sd)
{
    std::mt19937_64 rng(sd);
    const ll        R = static_cast<ll>(1) << 24;
    const size_t    n0 = 70000;
    std::vector<P>  pts = gen(rng, n0, R);
    {
        const std::vector<qdt_vec2d> d = to_d(pts);
        qdt_t* q = qdt_create();
        exec_fn(par_exec)(q, d.data(), n0);
        const InsFn ins = ins_fn(par_ins);
        std::string msg;
        for (int i = 0; i < 300 && msg.empty(); ++i)
        {
            P p;
            switch (rng() % 3)
            {
            case 0: p = pts[rng() % pts.size()]; break; /* duplicate */
            case 1:
            {
                const P a = pts[rng() % pts.size()], b = pts[rng() % pts.size()];
                p = P{ (a.x + b.x) / 2, (a.y + b.y) / 2 };
                if ((a.x + b.x) % 2 || (a.y + b.y) % 2) p = gen(rng, 1, R)[0];
                break;
            }
            default: p = gen(rng, 1, R)[0];
            }
            const size_t ret = ins(q, qdt_vec2d{ static_cast<double>(p.x), static_cast<double>(p.y) });
            pts.push_back(p);
            if (ret != pts.size())
            {
                msg = "insert returned the wrong count";
            }
        }
        if (msg.empty())
        {
            msg = verify(q, pts, ConSet(), false);
        }

        const std::vector<size_t> rep = lowest_dup(pts);
        std::vector<Edge>         acc;
        for (int i = 0; i < 14 && msg.empty(); ++i)
        {
            /* every other constraint starts at a point inserted after the mesh was built */
            const size_t a = (i % 2) ? pts.size() - 1 - rng() % 300 : rng() % pts.size();
            const size_t b = (i % 3 == 0) ? pts.size() - 1 - rng() % 300 : rng() % pts.size();
            const bool   ok = constraint_ok(pts, rep, acc, a, b, true);
            const TriList before = tri_sorted(tris_of(q));
            const size_t  rv = qdt_insert_constraint(q, a, b);
            if (rv != (ok ? 1u : 0u))
            {
                msg = fmt("%s constraint (%zu, %zu)", rv ? "unexpected accept of" : "wrongly rejected", a, b);
            }
            else if (!rv && tri_sorted(tris_of(q)) != before)
            {
                msg = "a rejected constraint changed the triangulation";
            }
            else if (ok)
            {
                acc.push_back(min_max(rep[a], rep[b]));
            }
        }
        if (msg.empty())
        {
            msg = verify(q, pts, ConSet(acc.begin(), acc.end()), false);
        }
        if (!msg.empty())
        {
            msg += fmt(" (%s batch call, %s insertions)", par_exec ? "parallel" : "serial", par_ins ? "parallel" : "serial");
        }
        r.finish(sd, msg);
        qdt_destroy(q);
    }
}

static void sc_parmesh(Scenario& r, const Config& c)
{
    const uint64_t seeds = c.quick ? 1 : 3;
    for (uint64_t s = 0; s < seeds; ++s)
    {
        for (int combo = 0; combo < 4; ++combo)
        {
            parmesh_case(r, (combo & 1) != 0, (combo & 2) != 0, seed_of(c, 81000 + 10 * s + static_cast<uint64_t>(combo)));
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Scenarios, part 3: large inputs, concurrent contexts                     */
/* ------------------------------------------------------------------------ */

static std::vector<P> big_points(size_t n, int dist, uint64_t sd)
{
    std::mt19937_64 rng(sd);
    std::vector<P>  p(n);
    for (size_t i = 0; i < n; ++i)
    {
        switch (dist)
        {
        case 0: /* uniform: general position for all practical purposes */
            p[i] = P{ static_cast<ll>(rng() % (1ULL << 28)), static_cast<ll>(rng() % (1ULL << 28)) };
            break;
        case 1: /* many duplicates and a cocircular grid */
            p[i] = P{ static_cast<ll>(rng() % 1500), static_cast<ll>(rng() % 1500) };
            break;
        default: /* flat band */
            p[i] = P{ static_cast<ll>(rng() % (1ULL << 28)), static_cast<ll>(rng() % 1024) };
            break;
        }
    }
    return p;
}

/* Inputs above the 524,288-point threshold, where the parallel call really
   splits the work: both results must be valid, and for inputs in general
   position they must be the same set of triangles. */
static void sc_big(Scenario& r, const Config& c)
{
    static const size_t sizes[3] = { 700000, 700000, 600000 };
    for (int dist = 0; dist < 3; ++dist)
    {
        const uint64_t        sd = seed_of(c, 4242 + static_cast<uint64_t>(dist));
        const std::vector<P>  pts = big_points(sizes[dist], dist, sd);
        const std::vector<qdt_vec2d> d = to_d(pts);
        qdt_t* a = qdt_create();
        qdt_t* b = qdt_create();
        qdt_execute(a, d.data(), d.size());
        qdt_execute_parallel(b, d.data(), d.size());
        std::string msg = verify(a, pts, ConSet(), false);
        if (!msg.empty())
        {
            msg = "serial: " + msg;
        }
        if (msg.empty())
        {
            msg = verify(b, pts, ConSet(), false);
            if (!msg.empty())
            {
                msg = "parallel: " + msg;
            }
        }
        if (msg.empty())
        {
            if (qdt_triangles_size(a) != qdt_triangles_size(b))
            {
                msg = fmt("serial and parallel triangle counts differ: %zu and %zu", qdt_triangles_size(a), qdt_triangles_size(b));
            }
            else if (dist != 1 && tri_sorted(tris_of(a)) != tri_sorted(tris_of(b)))
            {
                msg = "serial and parallel triangle sets differ";
            }
        }
        if (!msg.empty())
        {
            msg += fmt(" (%zu points, distribution %d)", pts.size(), dist);
        }
        std::fprintf(stderr, "  big: %zu points, distribution %d: %zu triangles\n", pts.size(), dist, qdt_triangles_size(a));
        r.finish(sd, msg);
        qdt_destroy(a);
        qdt_destroy(b);
    }
    r.note = "Distributions: uniform, a 1500 x 1500 grid with many duplicates (counts compared, not sets: cocircular ties), a flat band.";
}

static void contexts_worker(Scenario& r, const Config& c, unsigned t, uint64_t cases)
{
    static const ll ranges[] = { 3, 6, 20, 100, 10000 };
    for (uint64_t s = 0; s < cases; ++s)
    {
        const uint64_t  sd = seed_of(c, 200000 + 7919ULL * t + s);
        std::mt19937_64 rng(sd);
        const ll        R = ranges[rng() % (sizeof(ranges) / sizeof(ranges[0]))];
        const size_t    n = 3 + rng() % 60;
        std::vector<P>  pts = gen(rng, n, R);
        const bool      par = (s & 1) != 0;
        const std::vector<qdt_vec2d> d = to_d(pts);
        qdt_t* q = qdt_create();
        exec_fn(par)(q, d.data(), n);
        std::string msg = verify(q, pts, ConSet(), true);
        const size_t k = rng() % 4;
        for (size_t i = 0; i < k && msg.empty(); ++i)
        {
            const P p = gen(rng, 1, R)[0];
            ins_fn(!par)(q, qdt_vec2d{ static_cast<double>(p.x), static_cast<double>(p.y) });
            pts.push_back(p);
            msg = verify(q, pts, ConSet(), true);
        }
        r.finish(sd, msg);
        qdt_destroy(q);
    }
}

/* README: "Different contexts share no state and can be used concurrently." */
static void sc_contexts(Scenario& r, const Config& c)
{
    unsigned hw = std::thread::hardware_concurrency();
    if (hw < 2) hw = 2;
    if (hw > 8) hw = 8;
    {
        std::vector<std::thread> th;
        const uint64_t           cases = c.iters / 10 + 1;
        for (unsigned t = 0; t < hw; ++t)
        {
            th.push_back(std::thread([&r, &c, t, cases]() { contexts_worker(r, c, t, cases); }));
        }
        for (size_t i = 0; i < th.size(); ++i)
        {
            th[i].join();
        }
    }
    r.note = fmt("%u threads, each with its own contexts (batch, parallel batch and insertions).", hw);

    if (c.quick)
    {
        return;
    }
    /* Two threads run qdt_execute_parallel at the same time on one shared (read-only)
       input of 600,000 points; each result must equal the serial one. */
    const uint64_t              sd = seed_of(c, 4300);
    const std::vector<P>        pts = big_points(600000, 0, sd);
    const std::vector<qdt_vec2d> d = to_d(pts);
    qdt_t* ref = qdt_create();
    qdt_execute(ref, d.data(), d.size());
    const TriList want = tri_sorted(tris_of(ref));
    qdt_destroy(ref);
    qdt_t* q[2] = { qdt_create(), qdt_create() };
    {
        std::vector<std::thread> th;
        for (int i = 0; i < 2; ++i)
        {
            th.push_back(std::thread([&d, &q, i]() { qdt_execute_parallel(q[i], d.data(), d.size()); }));
        }
        for (size_t i = 0; i < th.size(); ++i)
        {
            th[i].join();
        }
    }
    for (int i = 0; i < 2; ++i)
    {
        r.finish(sd, tri_sorted(tris_of(q[i])) == want ? "" : "a context built while another one was running differs from the serial result");
        qdt_destroy(q[i]);
    }
    r.note += " Plus two parallel batch calls at the same time on 600,000 points.";
}

/* ------------------------------------------------------------------------ */
/* Scenario table                                                           */
/* ------------------------------------------------------------------------ */

struct Entry
{
    const char* name;
    const char* title;
    const char* section;
    void (*fn)(Scenario&, const Config&);
    bool heavy; /* needs inputs with hundreds of thousands of points; not run with --quick */
};

static const Entry kEntries[] = {
    { "api", "API contract: NULL context, out-of-range ids, query functions, duplicates, degenerate contexts, constraint end points, replacement by the next batch call",
      "[Functions](#functions), [Results](#results)", sc_api, false },
    { "examples", "The README examples (quadrilateral with a constraint, six-point program); the input array is copied and not modified",
      "[Example](#example), [Constrained triangulation](#constrained-triangulation)", sc_examples, false },
    { "guards", "Invalid input: NaN, infinity, NULL array, more than 100,000,000 points", "[Results](#results)", sc_guards, false },
    { "structured", "Degenerate and structured inputs: grids, collinear points, duplicates, fewer than three points, coordinates up to 2^28",
      "[Results](#results)", sc_structured, false },
    { "batch", "Random small point sets with many duplicates and cocircular points, serial and parallel batch call", "[Features](#features), [Results](#results)",
      sc_batch, false },
    { "magnitudes", "One point set that mixes magnitudes from 1 to 2^28, including clusters far from the origin", "[Results](#results)", sc_magnitudes, false },
    { "scale", "Scaling the whole input by 2^-1000 ... 2^900 must not change the result", "[Results](#results)", sc_scale, false },
    { "circle", "Nearly cocircular points, with insertions on top", "[Features](#features)", sc_circle, false },
    { "insert", "Insertion into small triangulations, also into a context that was never filled: duplicates, points on edges, serial and parallel",
      "[Batch and incremental usage](#batch-and-incremental-usage)", sc_insert, false },
    { "insertlarge", "Insertion into triangulations of a few thousand points: duplicates, points on edges, points outside the hull",
      "[Batch and incremental usage](#batch-and-incremental-usage)", sc_insertlarge, false },
    { "parmesh", "First insertion builds the internal mesh in parallel (70,000 points, at least 50,000 triangles), then insertions and constraints",
      "[Parallel mode](#parallel-mode)", sc_parmesh, false },
    { "cdt", "Constraints against an independent oracle for the rules: accepted, rejected, unchanged after a rejection; single calls and batch",
      "[Constrained triangulation](#constrained-triangulation)", sc_cdt, false },
    { "cdtins", "The same with insertions in between: constraints are kept, and split by points inserted on them",
      "[Constrained triangulation](#constrained-triangulation)", sc_cdtins, false },
    { "cdtequiv", "A batch call equals the single calls in the same order; in general position the order of the constraints does not matter",
      "[Constrained triangulation](#constrained-triangulation)", sc_cdtequiv, false },
    { "big", "600,000 to 700,000 points: serial and parallel batch call give valid results and the same triangles", "[Parallel mode](#parallel-mode)", sc_big, true },
    { "contexts", "Independent contexts used from several threads at the same time", "[Thread safety](#thread-safety)", sc_contexts, false },
};

static const size_t kEntryCount = sizeof(kEntries) / sizeof(kEntries[0]);

/* ------------------------------------------------------------------------ */
/* System information (same approach as test/qdt_bench.cpp)                 */
/* ------------------------------------------------------------------------ */

#if !defined(_WIN32)

/* Reads "key: value ..." from a /proc style file; returns the text after the colon. */
static bool proc_value(const char* path, const char* key, std::string& out)
{
    FILE* f = std::fopen(path, "r");
    if (!f)
    {
        return false;
    }
    char         line[512];
    const size_t klen = std::strlen(key);
    bool         ok = false;
    while (std::fgets(line, sizeof(line), f))
    {
        if (std::strncmp(line, key, klen) != 0)
        {
            continue;
        }
        const char* p = line + klen;
        while (*p == ' ' || *p == '\t')
        {
            ++p;
        }
        if (*p != ':')
        {
            continue;
        }
        ++p;
        out = p;
        trim_str(out);
        ok = true;
        break;
    }
    std::fclose(f);
    return ok;
}

#else

/* Reads a REG_SZ value from the 64-bit registry view. */
static bool reg_string(HKEY root, const char* subkey, const char* name, std::string& out)
{
    HKEY key = NULL;
    if (RegOpenKeyExA(root, subkey, 0, KEY_READ | KEY_WOW64_64KEY, &key) != ERROR_SUCCESS)
    {
        return false;
    }
    char buf[512];
    std::memset(buf, 0, sizeof(buf));
    DWORD size = static_cast<DWORD>(sizeof(buf) - 1);
    DWORD type = 0;
    const LONG rc = RegQueryValueExA(key, name, NULL, &type, reinterpret_cast<LPBYTE>(buf), &size);
    RegCloseKey(key);
    if (rc != ERROR_SUCCESS || type != REG_SZ)
    {
        return false;
    }
    out = buf;
    trim_str(out);
    return true;
}

#endif

static std::string cpu_name()
{
    std::string v;
#if defined(_WIN32)
    if (reg_string(HKEY_LOCAL_MACHINE, "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0", "ProcessorNameString", v) && !v.empty())
    {
        return v;
    }
#else
    if (proc_value("/proc/cpuinfo", "model name", v) && !v.empty())
    {
        return v;
    }
    if (proc_value("/proc/cpuinfo", "Model", v) && !v.empty())
    {
        return v;
    }
#endif
    return "unknown";
}

/* Installed physical memory in bytes, 0 if unknown. */
static unsigned long long total_memory_bytes()
{
#if defined(_WIN32)
    MEMORYSTATUSEX ms;
    std::memset(&ms, 0, sizeof(ms));
    ms.dwLength = sizeof(ms);
    if (GlobalMemoryStatusEx(&ms))
    {
        return static_cast<unsigned long long>(ms.ullTotalPhys);
    }
    return 0;
#else
    std::string v;
    if (proc_value("/proc/meminfo", "MemTotal", v))
    {
        return std::strtoull(v.c_str(), NULL, 10) * 1024ULL; /* the value is in kB */
    }
    return 0;
#endif
}

static std::string os_name()
{
#if defined(_WIN32)
    typedef LONG(WINAPI* RtlGetVersionFn)(PRTL_OSVERSIONINFOW);
    RTL_OSVERSIONINFOW vi;
    std::memset(&vi, 0, sizeof(vi));
    vi.dwOSVersionInfoSize = sizeof(vi);
    bool          have = false;
    const HMODULE nt = GetModuleHandleW(L"ntdll.dll");
    if (nt)
    {
        RtlGetVersionFn fn = reinterpret_cast<RtlGetVersionFn>(reinterpret_cast<void*>(GetProcAddress(nt, "RtlGetVersion")));
        if (fn && fn(&vi) == 0)
        {
            have = true;
        }
    }
    std::string name = "Windows";
    if (have)
    {
        if (vi.dwMajorVersion == 10)
        {
            name = (vi.dwBuildNumber >= 22000) ? "Windows 11" : "Windows 10";
        }
        else
        {
            name = fmt("Windows %lu.%lu", static_cast<unsigned long>(vi.dwMajorVersion), static_cast<unsigned long>(vi.dwMinorVersion));
        }
    }
    std::string display;
    if (reg_string(HKEY_LOCAL_MACHINE, "SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion", "DisplayVersion", display) && !display.empty())
    {
        name += " " + display;
    }
    if (have)
    {
        name += fmt(" (build %lu)", static_cast<unsigned long>(vi.dwBuildNumber));
    }
    SYSTEM_INFO si;
    std::memset(&si, 0, sizeof(si));
    GetNativeSystemInfo(&si);
    switch (si.wProcessorArchitecture)
    {
    case PROCESSOR_ARCHITECTURE_AMD64: name += ", x64"; break;
    case PROCESSOR_ARCHITECTURE_ARM64: name += ", ARM64"; break;
    case PROCESSOR_ARCHITECTURE_INTEL: name += ", x86"; break;
    default: break;
    }
    return name;
#elif defined(QDT_VERIFY_LINUX)
    struct utsname u;
    if (uname(&u) == 0)
    {
        return std::string(u.sysname) + " " + u.release + " " + u.machine;
    }
    return "unknown";
#else
    return "unknown";
#endif
}

static std::string compiler_string()
{
#if defined(__clang__)
    return std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("gcc ") + __VERSION__;
#elif defined(_MSC_VER)
#if defined(_MSC_FULL_VER)
    return fmt("MSVC %d.%d.%d", static_cast<int>(_MSC_FULL_VER / 10000000), static_cast<int>((_MSC_FULL_VER / 100000) % 100),
        static_cast<int>(_MSC_FULL_VER % 100000));
#else
    return fmt("MSVC (_MSC_VER %d)", static_cast<int>(_MSC_VER));
#endif
#else
    return "unknown";
#endif
}

/* ------------------------------------------------------------------------ */
/* Report                                                                   */
/* ------------------------------------------------------------------------ */

static const char kBegin[] = "<!-- QDT-VERIFY-BEGIN -->";
static const char kEnd[] = "<!-- QDT-VERIFY-END -->";

static std::string fmt_s(double ms)
{
    return fmt("%.2f", ms / 1000.0);
}

static std::string build_report(const Config& cfg, const std::vector<std::string>& notes, const std::deque<Scenario>& sc, double total_ms)
{
    unsigned long long cases = 0, failures = 0;
    size_t             skipped = 0;
    for (size_t i = 0; i < sc.size(); ++i)
    {
        cases += sc[i].cases;
        failures += sc[i].failures;
        if (sc[i].skipped) ++skipped;
    }
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 1;
    const unsigned long long mem = total_memory_bytes();
    char                     when[64];
    std::time_t              now = std::time(NULL);
    std::strftime(when, sizeof(when), "%Y-%m-%d %H:%M UTC", std::gmtime(&now));

    std::string o;
    o += kBegin;
    o += "\n## Correctness tests\n\n";
    o += "The results below were produced by `test/qdt_verify.cpp`, which uses only the public C API. It runs "
        "randomized and structured test cases and checks every triangulation with exact integer arithmetic, "
        "independent of the predicates inside the library: the inputs have integer coordinates with an absolute "
        "value of at most 2^28, which `double` represents exactly. A case is one input together with the "
        "operations applied to it (for example a point set, then insertions and constraints). The context is "
        "verified at the end of the case, and in the incremental scenarios also in between.\n\n";

    o += "### Test environment\n\n";
    o += std::string("- Date: ") + when + "\n";
    o += "- CPU: " + cpu_name() + "\n";
    o += fmt("- Hardware threads: %u\n", hw);
    if (mem != 0)
    {
        o += fmt("- Memory: %.1f GiB\n", static_cast<double>(mem) / (1024.0 * 1024.0 * 1024.0));
    }
    o += "- OS: " + os_name() + "\n";
    o += "- Compiler (test program): " + compiler_string() + fmt(", %zu-bit; exact arithmetic: ", sizeof(void*) * 8) + arith_name() + "\n";
    o += fmt("- Parameters: iterations = %llu, seed = %llu%s\n", static_cast<unsigned long long>(cfg.iters),
        static_cast<unsigned long long>(cfg.seed), cfg.quick ? ", quick mode" : "");
    for (size_t i = 0; i < notes.size(); ++i)
    {
        o += "- " + notes[i] + "\n";
    }
    o += "\n### Results\n\n";
    if (failures == 0)
    {
        o += fmt("**All %zu scenarios passed:** %s cases, 0 failures, %s s.", sc.size() - skipped, with_commas(cases).c_str(), fmt_s(total_ms).c_str());
        if (skipped != 0)
        {
            o += fmt(" %zu scenario(s) not run.", skipped);
        }
        if (cfg.partial)
        {
            o += " This was a partial run.";
        }
    }
    else
    {
        o += fmt("**%s of %s cases failed.** See the failures below.", with_commas(failures).c_str(), with_commas(cases).c_str());
    }
    o += "\n\n| Scenario | What is checked | README | Cases | Failures | Time (s) |\n";
    o += "| --- | --- | --- | ---: | ---: | ---: |\n";
    for (size_t i = 0; i < sc.size(); ++i)
    {
        const Scenario& s = sc[i];
        o += "| `" + s.name + "` | " + s.title + " | " + s.section + " | ";
        if (s.skipped)
        {
            o += "not run | - | - |\n";
        }
        else
        {
            o += with_commas(s.cases) + " | " + with_commas(s.failures) + " | " + fmt_s(s.ms) + " |\n";
        }
    }
    o += "\n";
    bool any_note = false;
    for (size_t i = 0; i < sc.size(); ++i)
    {
        if (!sc[i].note.empty())
        {
            if (!any_note) o += "Notes:\n\n";
            any_note = true;
            o += "- `" + sc[i].name + "`: " + sc[i].note + "\n";
        }
    }
    if (any_note) o += "\n";
    if (failures != 0)
    {
        o += "### Failures\n\n";
        for (size_t i = 0; i < sc.size(); ++i)
        {
            for (size_t k = 0; k < sc[i].samples.size(); ++k)
            {
                o += "- `" + sc[i].name + "`: " + sc[i].samples[k] + "\n";
            }
        }
        o += "\n";
    }

    o += "### What is checked in every case\n\n";
    o += "- Every triangle has three distinct vertices in strictly counter-clockwise order and refers only to the "
        "lowest index of a group of duplicate points.\n";
    o += "- Every directed edge occurs once, so triangles do not overlap; the areas of the triangles add up to the "
        "area of the convex hull of the points; no hull vertex lies to the right of a boundary edge.\n";
    o += "- Every interior edge is locally Delaunay, except for constrained edges. Without constraints, small "
        "inputs are also checked against the empty-circle property for every point.\n";
    o += "- Every accepted constraint is an edge of the triangulation, and a rejected constraint leaves the "
        "triangulation unchanged.\n";
    o += "- `qdt_points_size`, `qdt_get_point` and `qdt_point_is_used` agree with the points that were given and "
        "with the triangles.\n\n";

    o += "### Not covered\n\n";
    o += "- The fallback of `qdt_execute_parallel` to the serial algorithm (exhausted edge storage, threads that "
        "cannot be started) and out-of-memory behavior: these paths cannot be triggered through the public API.\n";
    o += "- Inputs of more than 700,000 points, and the speed and memory figures (see [Performance](#performance)).\n";
    o += "- Coordinates outside the range of the exact check: integers up to 2^28, scaled as a whole by powers of "
        "two. Subnormal coordinates are documented as unsupported.\n";
    o += "- Compiling `qdt.h` as C, and every compiler and platform other than the one listed above.\n";
    o += "- Using one context from several threads, which is not supported.\n\n";

    o += "### Reproducing the tests\n\n";
    o += "`test/qdt_verify.cpp` uses only the public API. `--readme README.md` "
        "replaces this section (the text between the two `QDT-VERIFY` marker comments, inserted before "
        "\"Current limitations\" the first time) and is refused when a check fails. The program exits with "
        "status 1 if a check fails. `--quick` runs a short smoke test and `--scenario LIST` selects scenarios; "
        "neither updates the README unless `--allow-partial` is given.\n\n";
    o += "```sh\n";
    o += "# Linux, with libqdt.a built as described in \"Building\"\n";
    o += "g++ -std=c++14 -O2 -DNDEBUG -pthread -Ilibqdt test/qdt_verify.cpp libqdt.a -lm -o qdt_verify\n";
    o += "./qdt_verify --readme README.md --note \"libqdt <commit>, g++ -O2, static library\"\n";
    o += "```\n\n";
    o += "```bat\n";
    o += ":: Windows, MSVC; qdt.lib built as a static library with /DQDT_DLL\n";
    o += "cl /std:c++14 /O2 /EHsc /DNDEBUG /DQDT_DLL /Ilibqdt test\\qdt_verify.cpp qdt.lib\n";
    o += "qdt_verify --readme README.md --note \"libqdt <commit>, MSVC /O2, static library\"\n";
    o += "```\n";
    o += kEnd;
    return o;
}

/* ------------------------------------------------------------------------ */
/* Files and README update                                                  */
/* ------------------------------------------------------------------------ */

static bool read_file(const std::string& path, std::string& out)
{
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f)
    {
        return false;
    }
    out.clear();
    char   buf[65536];
    size_t got;
    while ((got = std::fread(buf, 1, sizeof(buf), f)) > 0)
    {
        out.append(buf, got);
    }
    const bool ok = std::ferror(f) == 0;
    std::fclose(f);
    return ok;
}

static bool write_file(const std::string& path, const std::string& text)
{
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f)
    {
        return false;
    }
    const bool wrote = std::fwrite(text.data(), 1, text.size(), f) == text.size();
    const bool closed = std::fclose(f) == 0;
    return wrote && closed;
}

/* Position of the line that starts with the given heading, or npos. */
static size_t find_heading(const std::string& text, const char* heading)
{
    const std::string h(heading);
    size_t            pos = 0;
    while ((pos = text.find(h, pos)) != std::string::npos)
    {
        if (pos == 0 || text[pos - 1] == '\n')
        {
            return pos;
        }
        pos += h.size();
    }
    return std::string::npos;
}

/* Replaces the section between the markers, or inserts it the first time. */
static bool update_readme(const std::string& path, const std::string& block, std::string& what)
{
    std::string text;
    if (!read_file(path, text))
    {
        what = "cannot read " + path;
        return false;
    }
    const bool  crlf = text.find("\r\n") != std::string::npos;
    std::string blk;
    for (size_t i = 0; i < block.size(); ++i)
    {
        if (block[i] == '\n' && crlf) blk += '\r';
        blk += block[i];
    }
    const std::string nl = crlf ? "\r\n" : "\n";
    const size_t      b = text.find(kBegin);
    const size_t      e = text.find(kEnd);
    if (b != std::string::npos && e != std::string::npos && e > b)
    {
        text.replace(b, e + std::strlen(kEnd) - b, blk);
        what = "replaced the existing section";
    }
    else if (b != std::string::npos || e != std::string::npos)
    {
        what = path + " contains only one of the two QDT-VERIFY markers (or they are in the wrong order); not modified";
        return false;
    }
    else
    {
        size_t pos = find_heading(text, "## Current limitations");
        if (pos == std::string::npos)
        {
            pos = find_heading(text, "## License");
        }
        if (pos != std::string::npos)
        {
            text.insert(pos, blk + nl + nl);
            what = "inserted a new section";
        }
        else
        {
            if (!text.empty() && text[text.size() - 1] != '\n') text += nl;
            text += nl + blk + nl;
            what = "appended a new section at the end";
        }
    }
    if (!write_file(path, text))
    {
        what = "cannot write " + path;
        return false;
    }
    return true;
}

/* ------------------------------------------------------------------------ */
/* Driver                                                                   */
/* ------------------------------------------------------------------------ */

static void usage(const char* argv0)
{
    std::fprintf(stderr,
        "usage: %s [options]\n"
        "  --readme FILE      update the QDT-VERIFY section of FILE (only if every check passed)\n"
        "  --output FILE      write the Markdown section to FILE (UTF-8, LF)\n"
        "  --note TEXT        add a line to the \"Test environment\" list (repeatable)\n"
        "  --quick            fewer cases, no inputs with hundreds of thousands of points\n"
        "  --iterations N     base number of random cases per scenario (default 3000)\n"
        "  --seed N           offset of all random seeds (default 0)\n"
        "  --scenario LIST    comma separated scenario names, or \"all\" (default)\n"
        "  --allow-partial    allow --readme together with --quick or --scenario\n"
        "  --list             list the scenarios\n"
        "  --help\n"
        "Without --readme and --output the section is printed to stdout; progress goes to stderr.\n"
        "Exit status: 0 all passed, 1 a check failed, 2 usage or file error.\n",
        argv0);
}

static bool parse_list(const std::string& list, std::set<std::string>& out)
{
    size_t pos = 0;
    while (pos <= list.size())
    {
        size_t comma = list.find(',', pos);
        if (comma == std::string::npos) comma = list.size();
        const std::string item = list.substr(pos, comma - pos);
        pos = comma + 1;
        if (item.empty()) continue;
        if (item == "all") continue;
        bool known = false;
        for (size_t i = 0; i < kEntryCount; ++i)
        {
            if (item == kEntries[i].name) known = true;
        }
        if (!known)
        {
            std::fprintf(stderr, "unknown scenario: %s\n", item.c_str());
            return false;
        }
        out.insert(item);
    }
    return true;
}

int main(int argc, char** argv)
{
    Config                   cfg;
    std::string              readme_path, output_path;
    std::vector<std::string> notes;
    std::set<std::string>    only;
    bool                     allow_partial = false, iters_set = false;

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--quick")
        {
            cfg.quick = true;
        }
        else if (a == "--allow-partial")
        {
            allow_partial = true;
        }
        else if (a == "--iterations" && i + 1 < argc)
        {
            cfg.iters = static_cast<uint64_t>(std::strtoull(argv[++i], NULL, 10));
            iters_set = true;
            if (cfg.iters == 0)
            {
                std::fprintf(stderr, "--iterations needs a positive number\n");
                return 2;
            }
        }
        else if (a == "--seed" && i + 1 < argc)
        {
            cfg.seed = static_cast<uint64_t>(std::strtoull(argv[++i], NULL, 10));
        }
        else if (a == "--scenario" && i + 1 < argc)
        {
            if (!parse_list(argv[++i], only)) return 2;
        }
        else if (a == "--readme" && i + 1 < argc)
        {
            readme_path = argv[++i];
        }
        else if (a == "--output" && i + 1 < argc)
        {
            output_path = argv[++i];
        }
        else if (a == "--note" && i + 1 < argc)
        {
            notes.push_back(argv[++i]);
        }
        else if (a == "--list")
        {
            for (size_t k = 0; k < kEntryCount; ++k)
            {
                std::printf("%-12s %s\n", kEntries[k].name, kEntries[k].title);
            }
            return 0;
        }
        else
        {
            usage(argv[0]);
            return a == "--help" ? 0 : 2;
        }
    }
    if (cfg.quick && !iters_set)
    {
        cfg.iters = 300;
    }
    cfg.partial = cfg.quick || (!only.empty() && only.size() < kEntryCount);

    std::deque<Scenario> sc;
    const double         t_all = now_ms();
    for (size_t k = 0; k < kEntryCount; ++k)
    {
        const Entry& e = kEntries[k];
        if (!only.empty() && !only.count(e.name)) continue;
        sc.emplace_back();
        Scenario& s = sc.back();
        s.name = e.name;
        s.title = e.title;
        s.section = e.section;
        if (cfg.quick && e.heavy)
        {
            s.skipped = true;
            s.skip_reason = "--quick";
            std::fprintf(stderr, "[%s] not run (--quick)\n", e.name);
            continue;
        }
        std::fprintf(stderr, "[%s] running ...\n", e.name);
        const double t0 = now_ms();
        e.fn(s, cfg);
        s.ms = now_ms() - t0;
        std::fprintf(stderr, "[%s] %llu cases, %llu failures, %s s\n", e.name, s.cases, s.failures, fmt_s(s.ms).c_str());
    }
    const double total_ms = now_ms() - t_all;

    unsigned long long failures = 0;
    for (size_t i = 0; i < sc.size(); ++i) failures += sc[i].failures;

    const std::string report = build_report(cfg, notes, sc, total_ms);
    int               status = failures ? 1 : 0;

    if (output_path.empty() && readme_path.empty())
    {
        std::fputs(report.c_str(), stdout);
        std::fputc('\n', stdout);
    }
    if (!output_path.empty())
    {
        if (!write_file(output_path, report + "\n"))
        {
            std::fprintf(stderr, "cannot write %s\n", output_path.c_str());
            return 2;
        }
        std::fprintf(stderr, "report written to %s\n", output_path.c_str());
    }
    if (!readme_path.empty())
    {
        if (failures != 0)
        {
            std::fprintf(stderr, "%s not modified: %llu check(s) failed\n", readme_path.c_str(), failures);
        }
        else if (cfg.partial && !allow_partial)
        {
            std::fprintf(stderr, "%s not modified: this was a partial run (--quick or --scenario); use --allow-partial to write it anyway\n",
                readme_path.c_str());
            status = 2;
        }
        else
        {
            std::string what;
            if (!update_readme(readme_path, report, what))
            {
                std::fprintf(stderr, "%s\n", what.c_str());
                return 2;
            }
            std::fprintf(stderr, "%s: %s\n", readme_path.c_str(), what.c_str());
        }
    }
    if (failures)
    {
        std::fprintf(stderr, "FAILED: %llu check(s) failed\n", failures);
    }
    else
    {
        std::fprintf(stderr, "all checks passed\n");
    }
    return status;
}