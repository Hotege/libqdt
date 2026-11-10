/*
 * qdt_bench.cpp - performance benchmark for libqdt.
 *
 * The program only uses the public C API in qdt.h and builds on Linux and
 * Windows (MSVC, MinGW-w64, GCC, Clang). It produces a report in
 * GitHub-flavoured Markdown that is meant to be pasted into README.md as it
 * is: the report starts with a "## Performance" heading, its sections are
 * "###" headings and the per-insert tables "####" headings, and nothing but
 * Markdown goes into it. Progress messages go to stderr. Optionally every
 * measurement is also written to a CSV file.
 *
 * Scenarios (select with --scenario, default: all):
 *
 *   batch       qdt_execute / qdt_execute_parallel on uniform random points,
 *               100 ... 2,624,409 points (the size used by the old README)
 *   dist        the same calls on different point distributions: uniform,
 *               presorted, gaussian, clustered, regular grid, jittered grid,
 *               points on a circle, flat band, 50% duplicates
 *   insert      qdt_insert_point / qdt_insert_point_parallel: cost of the
 *               first insertion (internal mesh build), later insertions, the
 *               break-even point against a full rebuild, per-insert detail
 *   constraint  constrained Delaunay: qdt_insert_constraint and
 *               qdt_insert_constraints_batch with short non-crossing
 *               constraints, long chords across the whole domain, and a
 *               closed polygon ring
 *   query       qdt_get_triangle throughput
 *   concurrent  independent contexts triangulating at the same time on
 *               1, 2, 4 ... hardware threads (scaling and thread safety)
 *   memory      measured memory use per point (Linux and Windows); every
 *               measurement runs in a process of its own
 *   cgal        optional comparison with CGAL (build with -DQDT_BENCH_CGAL)
 *
 * Every scenario also performs cheap sanity checks (serial and parallel
 * results agree, insertions return the expected point count, accepted
 * constraints do not change the triangle count). If one fails, a warning block
 * is added to the end of the report and the exit code is 1.
 *
 * Build on Linux / macOS (the library must have been built with
 * -std=c++14 or newer and with optimizations, see the repository README):
 *
 *   g++ -std=c++14 -O2 -DNDEBUG -pthread -I../libqdt qdt_bench.cpp \
 *       <path-to>/libqdt.a -lm -o qdt_bench
 *
 * Build on Windows with MSVC (x64 Native Tools Command Prompt). qdt.h marks
 * the API __declspec(dllimport) unless QDT_DLL is defined (see "Windows
 * notes" in the README), so the command line depends on how libqdt was built.
 * Both variants need /O2 (or /Ox) and no /fp:fast:
 *
 *   static library (qdt.lib compiled with /DQDT_DLL, as the README says):
 *     cl /std:c++14 /O2 /EHsc /DNDEBUG /DQDT_DLL /I..\libqdt qdt_bench.cpp qdt.lib
 *
 *   DLL (qdt.lib is the import library; qdt.dll must be next to the .exe):
 *     cl /std:c++14 /O2 /EHsc /DNDEBUG /I..\libqdt qdt_bench.cpp qdt.lib
 *
 * Build on Windows with MinGW-w64 (posix thread model; static library built
 * with -DQDT_DLL):
 *
 *   x86_64-w64-mingw32-g++ -std=c++14 -O2 -DNDEBUG -DQDT_DLL -pthread \
 *       -I../libqdt qdt_bench.cpp <path-to>/libqdt.a -lpsapi -static \
 *       -o qdt_bench.exe
 *
 * Build with the CGAL comparison (CGAL 5.x or newer, headers and GMP/MPFR
 * installed; link GMP/MPFR only if your CGAL installation uses them):
 *
 *   g++ -std=c++17 -O2 -DNDEBUG -DQDT_BENCH_CGAL -pthread -I../libqdt \
 *       qdt_bench.cpp <path-to>/libqdt.a -lgmp -lmpfr -lm -o qdt_bench
 *
 * Usage:
 *   qdt_bench [options] --output bench.md
 *     --output FILE        write the Markdown report to FILE (UTF-8, LF line
 *                          ends) instead of stdout. Prefer this to shell
 *                          redirection: Windows PowerShell 5.1 writes
 *                          redirected output as UTF-16, which breaks Markdown
 *     --note TEXT          add a line to the "Test environment" list, for
 *                          example the libqdt commit or how the library was
 *                          built; can be given several times
 *     --quick              small sizes only (a short smoke run)
 *     --repeat N           repetitions per measurement (default 5; very large
 *                          inputs and the insert/constraint scenarios use at
 *                          most 3)
 *     --seed N             seed of the pseudo random generator
 *     --scenario LIST      comma separated scenarios, or "all"
 *     --max-points N       skip every size above N
 *     --csv FILE           also write all measurements to FILE
 *     --help
 *
 * For numbers that are worth publishing: use an optimized build of the
 * library, close other programs, keep the machine on mains power with a
 * high-performance power plan / fixed CPU governor, and run the full benchmark
 * once without --quick. A full run takes several minutes. The first-insertion
 * time is the least stable number (it is dominated by fresh memory being
 * mapped); repeat the run once or twice before quoting it.
 *
 * All timings are wall-clock milliseconds measured with std::chrono::
 * steady_clock around the library calls only; point generation, context
 * creation and destruction are not timed. Tables show medians.
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
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#if defined(__linux__)
#include <sys/types.h>
#include <sys/utsname.h>
#include <sys/wait.h>
#include <unistd.h>
#if defined(__GLIBC__)
#include <malloc.h>
#endif
#define QDT_BENCH_LINUX 1
#endif

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#if defined(_MSC_VER)
#pragma comment(lib, "psapi.lib")
#endif
#endif

#if defined(QDT_BENCH_LINUX) || defined(_WIN32)
#define QDT_BENCH_HAVE_MEMORY 1
#endif

#if defined(QDT_BENCH_CGAL)
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/version.h>
#endif

 /* ------------------------------------------------------------------------ */
 /* Basics: timing, formatting, statistics, random numbers                   */
 /* ------------------------------------------------------------------------ */

typedef std::chrono::steady_clock Clock;
typedef std::vector<qdt_vec2d>    PointList;

static const double kExtent = 1000.0;  /* points live in [0, kExtent)^2 */
static const double kPi = 3.14159265358979323846;
static const size_t kReadmeN = 2624409; /* point count of the old README benchmark */

static double now_ms()
{
    return std::chrono::duration<double, std::milli>(Clock::now().time_since_epoch()).count();
}

static std::string fmt(const char* f, ...)
{
    char    buf[512];
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

static std::string fmt_ms(double ms)
{
    if (ms < 1.0) return fmt("%.4f", ms);
    if (ms < 10.0) return fmt("%.3f", ms);
    if (ms < 100.0) return fmt("%.2f", ms);
    return fmt("%.1f", ms);
}

static std::string fmt_us(double us)
{
    if (us < 10.0) return fmt("%.2f", us);
    if (us < 100.0) return fmt("%.1f", us);
    return fmt("%.0f", us);
}

static std::string fmt_x(double x)
{
    return fmt("%.2fx", x);
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

struct Stats
{
    double min, median, mean, p95, max;
};

static Stats summarize(std::vector<double> v)
{
    Stats s;
    s.min = s.median = s.mean = s.p95 = s.max = 0.0;
    if (v.empty())
    {
        return s;
    }
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    s.min = v[0];
    s.max = v[n - 1];
    s.median = (n % 2 == 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sum += v[i];
    }
    s.mean = sum / static_cast<double>(n);
    size_t rank = static_cast<size_t>(std::ceil(0.95 * static_cast<double>(n)));
    if (rank < 1)
    {
        rank = 1;
    }
    s.p95 = v[rank - 1];
    return s;
}

/* splitmix64: small, fast and identical on every platform. */
struct Rng
{
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}

    uint64_t next()
    {
        uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    /* Uniform double in [0, 1). */
    double uni() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }

    size_t index(size_t n) { return static_cast<size_t>(next() % static_cast<uint64_t>(n)); }

    /* Standard normal (Box-Muller, the second value is discarded). */
    double gauss()
    {
        double u1 = uni();
        if (u1 < 1e-300)
        {
            u1 = 1e-300;
        }
        const double u2 = uni();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
    }
};

template <class T>
static void shuffle_vec(std::vector<T>& v, Rng& rng)
{
    for (size_t i = v.size(); i > 1; --i)
    {
        std::swap(v[i - 1], v[rng.index(i)]);
    }
}

static qdt_vec2d P(double x, double y)
{
    qdt_vec2d p;
    p.x = x;
    p.y = y;
    return p;
}

/* The C API takes a non-const pointer but never modifies the points. */
static qdt_vec2d* pdata(const PointList& v)
{
    static qdt_vec2d dummy = { 0.0, 0.0 };
    return v.empty() ? &dummy : const_cast<qdt_vec2d*>(&v[0]);
}

/* ------------------------------------------------------------------------ */
/* Report helpers                                                           */
/* ------------------------------------------------------------------------ */

static std::vector<std::string> g_warnings;
static volatile uint64_t        g_sink = 0; /* keeps read loops from being optimized away */

static void warn(const std::string& msg)
{
    g_warnings.push_back(msg);
    std::fprintf(stderr, "  [bench] WARNING: %s\n", msg.c_str());
    std::fflush(stderr);
}

static void progress(const std::string& msg)
{
    std::fprintf(stderr, "  [bench] %s\n", msg.c_str());
    std::fflush(stderr);
}

/* "###" heading plus an optional one-paragraph description. */
static void section(const char* title, const std::string& note)
{
    std::printf("### %s\n\n", title);
    if (!note.empty())
    {
        std::printf("%s\n\n", note.c_str());
    }
    std::fflush(stdout);
}

static void subsection(const std::string& title)
{
    std::printf("#### %s\n\n", title.c_str());
    std::fflush(stdout);
}

struct Row
{
    std::vector<std::string> cells;
    Row& operator()(const std::string& s)
    {
        cells.push_back(s);
        return *this;
    }
};

/* Markdown table. The first column and columns added with lcol() are left aligned, the others right aligned. */
class Table
{
public:
    /* Right aligned column (numbers). */
    Table& col(const std::string& h)
    {
        headers_.push_back(h);
        left_.push_back(headers_.size() == 1);
        return *this;
    }

    /* Left aligned column (text). */
    Table& lcol(const std::string& h)
    {
        headers_.push_back(h);
        left_.push_back(true);
        return *this;
    }

    void add(const Row& r) { rows_.push_back(r.cells); }

    bool empty() const { return rows_.empty(); }

    void print() const
    {
        if (rows_.empty())
        {
            std::printf("(no data)\n\n");
            return;
        }
        const size_t        cols = headers_.size();
        std::vector<size_t> w(cols, 3);
        for (size_t c = 0; c < cols; ++c)
        {
            w[c] = std::max(w[c], headers_[c].size());
        }
        for (size_t r = 0; r < rows_.size(); ++r)
        {
            for (size_t c = 0; c < cols && c < rows_[r].size(); ++c)
            {
                w[c] = std::max(w[c], rows_[r][c].size());
            }
        }

        print_row(headers_, w);
        std::printf("|");
        for (size_t c = 0; c < cols; ++c)
        {
            if (left_[c])
            {
                std::printf(" %s |", std::string(w[c], '-').c_str());
            }
            else
            {
                std::printf(" %s: |", std::string(w[c] - 1, '-').c_str());
            }
        }
        std::printf("\n");
        for (size_t r = 0; r < rows_.size(); ++r)
        {
            print_row(rows_[r], w);
        }
        std::printf("\n");
        std::fflush(stdout);
    }

private:
    void print_row(const std::vector<std::string>& r, const std::vector<size_t>& w) const
    {
        std::printf("|");
        for (size_t c = 0; c < w.size(); ++c)
        {
            const std::string s = (c < r.size()) ? r[c] : std::string();
            const size_t      pad = (w[c] > s.size()) ? (w[c] - s.size()) : 0;
            if (left_[c])
            {
                std::printf(" %s%*s |", s.c_str(), static_cast<int>(pad), "");
            }
            else
            {
                std::printf(" %*s%s |", static_cast<int>(pad), "", s.c_str());
            }
        }
        std::printf("\n");
    }

    std::vector<bool>                     left_;
    std::vector<std::string>              headers_;
    std::vector<std::vector<std::string>> rows_;
};

struct CsvRow
{
    std::string scenario;
    std::string name;
    size_t      n;
    std::string metric;
    double      value;
};

static std::vector<CsvRow> g_csv;

static void rec(const char* scenario, const std::string& name, size_t n, const char* metric, double value)
{
    CsvRow r;
    r.scenario = scenario;
    r.name = name;
    r.n = n;
    r.metric = metric;
    r.value = value;
    g_csv.push_back(r);
}

static std::string csv_quote(const std::string& s)
{
    std::string out = "\"";
    for (size_t i = 0; i < s.size(); ++i)
    {
        if (s[i] == '"')
        {
            out += "\"\"";
        }
        else
        {
            out += s[i];
        }
    }
    out += "\"";
    return out;
}

static bool write_csv(const std::string& path)
{
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f)
    {
        return false;
    }
    std::fprintf(f, "scenario,case,points,metric,value\n");
    for (size_t i = 0; i < g_csv.size(); ++i)
    {
        const CsvRow& r = g_csv[i];
        std::fprintf(f, "%s,%s,%zu,%s,%.9g\n", csv_quote(r.scenario).c_str(), csv_quote(r.name).c_str(), r.n,
            csv_quote(r.metric).c_str(), r.value);
    }
    return std::fclose(f) == 0;
}

/* ------------------------------------------------------------------------ */
/* Configuration                                                            */
/* ------------------------------------------------------------------------ */

struct Config
{
    bool                     quick;
    size_t                   reps;
    uint64_t                 seed;
    size_t                   max_points;
    std::string              csv_path;
    std::string              output_path;
    std::vector<std::string> notes;
    std::set<std::string>    scenarios;

    std::vector<size_t> batch_sizes;
    std::vector<size_t> insert_sizes;
    std::vector<size_t> constraint_sizes;
    std::vector<size_t> memory_sizes;
    std::vector<size_t> cgal_sizes;
    size_t              dist_n;
    size_t              query_n;
    size_t              concurrent_n;
    size_t              insert_count;
    size_t              short_constraints;
    size_t              long_chords;
    size_t              ring_points;

    Config()
        : quick(false), reps(5), seed(20260607ULL), max_points(0), dist_n(0), query_n(0), concurrent_n(0),
        insert_count(0), short_constraints(0), long_chords(0), ring_points(0)
    {
    }
};

static void limit_sizes(std::vector<size_t>& v, size_t max_points)
{
    if (max_points == 0)
    {
        return;
    }
    std::vector<size_t> out;
    for (size_t i = 0; i < v.size(); ++i)
    {
        if (v[i] <= max_points)
        {
            out.push_back(v[i]);
        }
    }
    v.swap(out);
}

static size_t cap_size(size_t v, size_t max_points)
{
    return (max_points != 0 && v > max_points) ? max_points : v;
}

static void build_plan(Config& c)
{
    if (c.quick)
    {
        const size_t b[] = { 100, 1000, 10000, 100000 };
        c.batch_sizes.assign(b, b + sizeof(b) / sizeof(b[0]));
        const size_t i[] = { 10000, 100000 };
        c.insert_sizes.assign(i, i + sizeof(i) / sizeof(i[0]));
        const size_t k[] = { 10000, 100000 };
        c.constraint_sizes.assign(k, k + sizeof(k) / sizeof(k[0]));
        const size_t m[] = { 100000 };
        c.memory_sizes.assign(m, m + sizeof(m) / sizeof(m[0]));
        const size_t g[] = { 10000, 100000 };
        c.cgal_sizes.assign(g, g + sizeof(g) / sizeof(g[0]));
        c.dist_n = 100000;
        c.query_n = 100000;
        c.concurrent_n = 20000;
        c.insert_count = 20;
        c.short_constraints = 5000;
        c.long_chords = 20;
        c.ring_points = 500;
    }
    else
    {
        const size_t b[] = { 100, 1000, 10000, 100000, 1000000, kReadmeN };
        c.batch_sizes.assign(b, b + sizeof(b) / sizeof(b[0]));
        const size_t i[] = { 10000, 100000, 1000000, kReadmeN };
        c.insert_sizes.assign(i, i + sizeof(i) / sizeof(i[0]));
        const size_t k[] = { 10000, 100000, 1000000 };
        c.constraint_sizes.assign(k, k + sizeof(k) / sizeof(k[0]));
        const size_t m[] = { 100000, 1000000, kReadmeN };
        c.memory_sizes.assign(m, m + sizeof(m) / sizeof(m[0]));
        const size_t g[] = { 10000, 100000, 1000000, kReadmeN };
        c.cgal_sizes.assign(g, g + sizeof(g) / sizeof(g[0]));
        c.dist_n = 1000000;
        c.query_n = 1000000;
        c.concurrent_n = 100000;
        c.insert_count = 50;
        c.short_constraints = 20000;
        c.long_chords = 100;
        c.ring_points = 2000;
    }

    limit_sizes(c.batch_sizes, c.max_points);
    limit_sizes(c.insert_sizes, c.max_points);
    limit_sizes(c.constraint_sizes, c.max_points);
    limit_sizes(c.memory_sizes, c.max_points);
    limit_sizes(c.cgal_sizes, c.max_points);
    c.dist_n = cap_size(c.dist_n, c.max_points);
    c.query_n = cap_size(c.query_n, c.max_points);
    c.concurrent_n = cap_size(c.concurrent_n, c.max_points);
}

/* Very large inputs are repeated fewer times to keep the total run time sane. */
static size_t reps_for(const Config& c, size_t n)
{
    return (n >= 2000000 && c.reps > 3) ? 3 : c.reps;
}

/* ------------------------------------------------------------------------ */
/* Point generators                                                         */
/* ------------------------------------------------------------------------ */

typedef PointList(*GenFn)(size_t, Rng&);

static PointList gen_uniform(size_t n, Rng& r)
{
    PointList pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        const double x = kExtent * r.uni();
        const double y = kExtent * r.uni();
        pts.push_back(P(x, y));
    }
    return pts;
}

static PointList gen_sorted(size_t n, Rng& r)
{
    PointList pts = gen_uniform(n, r);
    std::sort(pts.begin(), pts.end(), [](const qdt_vec2d& a, const qdt_vec2d& b) {
        if (a.x != b.x) return a.x < b.x;
        return a.y < b.y;
        });
    return pts;
}

static PointList gen_gaussian(size_t n, Rng& r)
{
    PointList pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        const double x = kExtent * 0.5 + kExtent * 0.1 * r.gauss();
        const double y = kExtent * 0.5 + kExtent * 0.1 * r.gauss();
        pts.push_back(P(x, y));
    }
    return pts;
}

static PointList gen_clustered(size_t n, Rng& r)
{
    const int K = 20;
    double    cx[K], cy[K], sg[K];
    for (int k = 0; k < K; ++k)
    {
        cx[k] = kExtent * (0.1 + 0.8 * r.uni());
        cy[k] = kExtent * (0.1 + 0.8 * r.uni());
        sg[k] = 5.0 + 35.0 * r.uni();
    }
    PointList pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        const size_t k = r.index(static_cast<size_t>(K));
        const double x = cx[k] + sg[k] * r.gauss();
        const double y = cy[k] + sg[k] * r.gauss();
        pts.push_back(P(x, y));
    }
    return pts;
}

/* Integer lattice: a huge number of cocircular quadruples. Shuffled. */
static PointList gen_grid(size_t n, Rng& r)
{
    const size_t w = static_cast<size_t>(std::ceil(std::sqrt(static_cast<double>(n))));
    PointList    pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        pts.push_back(P(static_cast<double>(i % w), static_cast<double>(i / w)));
    }
    shuffle_vec(pts, r);
    return pts;
}

/* Lattice with small random offsets: regular, mesh-like, few exact degeneracies. */
static PointList gen_jitter(size_t n, Rng& r)
{
    const size_t w = static_cast<size_t>(std::ceil(std::sqrt(static_cast<double>(n))));
    PointList    pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        const double jx = (r.uni() - 0.5) * 0.6;
        const double jy = (r.uni() - 0.5) * 0.6;
        pts.push_back(P(static_cast<double>(i % w) + jx, static_cast<double>(i / w) + jy));
    }
    shuffle_vec(pts, r);
    return pts;
}

/* Every point on one circle: all points are on the convex hull and nearly cocircular. */
static PointList gen_circle(size_t n, Rng& r)
{
    PointList pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        const double a = 2.0 * kPi * r.uni();
        pts.push_back(P(kExtent * 0.5 + kExtent * 0.5 * std::cos(a), kExtent * 0.5 + kExtent * 0.5 * std::sin(a)));
    }
    return pts;
}

/* Very thin strip (aspect ratio 1e6 : 1): nearly collinear, stresses the exact predicates. */
static PointList gen_band(size_t n, Rng& r)
{
    PointList pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        const double x = kExtent * r.uni();
        const double y = 1e-3 * r.uni();
        pts.push_back(P(x, y));
    }
    return pts;
}

/* Half of the points are exact copies of the other half. */
static PointList gen_dups(size_t n, Rng& r)
{
    const size_t m = (n + 1) / 2;
    PointList    pts = gen_uniform(m, r);
    pts.reserve(n);
    while (pts.size() < n)
    {
        pts.push_back(pts[r.index(m)]);
    }
    shuffle_vec(pts, r);
    return pts;
}

struct Dist
{
    const char* name;
    GenFn       gen;
};

static const Dist kDists[] = {
    { "uniform random", gen_uniform },
    { "uniform, presorted by x", gen_sorted },
    { "gaussian", gen_gaussian },
    { "clustered (20 gaussians)", gen_clustered },
    { "regular grid (shuffled)", gen_grid },
    { "jittered grid (shuffled)", gen_jitter },
    { "all points on a circle", gen_circle },
    { "flat band (aspect 1e6:1)", gen_band },
    { "50% duplicate points", gen_dups },
};
static const size_t kDistCount = sizeof(kDists) / sizeof(kDists[0]);

/* ------------------------------------------------------------------------ */
/* Measuring batch triangulation                                            */
/* ------------------------------------------------------------------------ */

/*
 * Runs `iters` consecutive batch calls on one fresh context and returns the
 * average time of one call in ms. Small inputs use several iterations so that
 * the timer resolution does not matter.
 */
static double time_execute(const PointList& pts, bool parallel, size_t iters, size_t* tris_out)
{
    qdt_t* q = qdt_create();
    qdt_vec2d* p = pdata(pts);
    const double t0 = now_ms();
    for (size_t i = 0; i < iters; ++i)
    {
        if (parallel)
        {
            qdt_execute_parallel(q, p, pts.size());
        }
        else
        {
            qdt_execute(q, p, pts.size());
        }
    }
    const double t1 = now_ms();
    if (tris_out)
    {
        *tris_out = qdt_triangles_size(q);
    }
    qdt_destroy(q);
    return (t1 - t0) / static_cast<double>(iters);
}

static size_t iters_for(size_t n)
{
    if (n >= 50000)
    {
        return 1;
    }
    size_t it = 200000 / (n > 0 ? n : 1);
    if (it < 1) it = 1;
    if (it > 2000) it = 2000;
    return it;
}

struct BatchResult
{
    Stats  ms;
    size_t tris;
};

static BatchResult measure_batch(const PointList& pts, bool parallel, size_t reps)
{
    const size_t n = pts.size();
    const size_t iters = iters_for(n);
    size_t       tris = 0;
    if (n <= 1000000)
    {
        time_execute(pts, parallel, iters, &tris); /* warm-up */
    }
    std::vector<double> t;
    for (size_t i = 0; i < reps; ++i)
    {
        t.push_back(time_execute(pts, parallel, iters, &tris));
    }
    BatchResult r;
    r.ms = summarize(t);
    r.tris = tris;
    return r;
}

static void bench_batch(const Config& c)
{
    section("Batch triangulation: scaling with the number of points",
        "Uniform random points; medians of " + std::to_string(c.reps) +
        " runs (3 for the largest sizes). `qdt_execute_parallel` uses several threads only from 200,000 "
        "distinct points on; below that it runs on the calling thread. Inputs below 50,000 points are "
        "timed as the average of several consecutive calls on one context.");

    Table t;
    t.col("Points").col("Triangles").col("Serial (ms)").col("Parallel (ms)").col("Speedup").col("Serial (Mpts/s)").col(
        "Parallel (Mpts/s)");

    for (size_t k = 0; k < c.batch_sizes.size(); ++k)
    {
        const size_t n = c.batch_sizes[k];
        progress("batch, " + with_commas(n) + " points");
        Rng       rng(c.seed + n);
        PointList pts = gen_uniform(n, rng);
        const size_t reps = reps_for(c, n);

        const BatchResult s = measure_batch(pts, false, reps);
        const BatchResult p = measure_batch(pts, true, reps);
        if (s.tris != p.tris)
        {
            warn(fmt("batch, %zu points: serial and parallel triangle counts differ (%zu vs %zu)", n, s.tris, p.tris));
        }

        const double sm = s.ms.median, pm = p.ms.median;
        t.add(Row()(with_commas(n))(with_commas(s.tris))(fmt_ms(sm))(fmt_ms(pm))(fmt_x(sm / pm))(
            fmt("%.2f", static_cast<double>(n) / sm / 1000.0))(fmt("%.2f", static_cast<double>(n) / pm / 1000.0)));

        rec("batch", "uniform", n, "triangles", static_cast<double>(s.tris));
        rec("batch", "uniform", n, "serial_ms_median", sm);
        rec("batch", "uniform", n, "serial_ms_min", s.ms.min);
        rec("batch", "uniform", n, "parallel_ms_median", pm);
        rec("batch", "uniform", n, "parallel_ms_min", p.ms.min);
    }
    t.print();
}

static void bench_dist(const Config& c)
{
    const size_t n = c.dist_n;
    section("Batch triangulation: point distributions",
        with_commas(n) + " points per case; medians of " + std::to_string(reps_for(c, n)) +
        " runs. Triangle counts differ between cases because the number of convex hull points and of "
        "duplicate points differs.");
    if (n == 0)
    {
        std::printf("(skipped: --max-points is too small)\n\n");
        return;
    }

    Table t;
    t.lcol("Distribution").col("Triangles").col("Serial (ms)").col("Parallel (ms)").col("Speedup").col("Serial (Mpts/s)");

    for (size_t d = 0; d < kDistCount; ++d)
    {
        progress(std::string("distribution: ") + kDists[d].name);
        Rng       rng(c.seed + 1000 + d);
        PointList pts = kDists[d].gen(n, rng);
        const size_t reps = reps_for(c, n);

        const BatchResult s = measure_batch(pts, false, reps);
        const BatchResult p = measure_batch(pts, true, reps);
        if (s.tris != p.tris)
        {
            warn(fmt("distribution '%s': serial and parallel triangle counts differ (%zu vs %zu)", kDists[d].name,
                s.tris, p.tris));
        }

        const double sm = s.ms.median, pm = p.ms.median;
        t.add(Row()(kDists[d].name)(with_commas(s.tris))(fmt_ms(sm))(fmt_ms(pm))(fmt_x(sm / pm))(
            fmt("%.2f", static_cast<double>(n) / sm / 1000.0)));

        rec("dist", kDists[d].name, n, "triangles", static_cast<double>(s.tris));
        rec("dist", kDists[d].name, n, "serial_ms_median", sm);
        rec("dist", kDists[d].name, n, "parallel_ms_median", pm);
    }
    t.print();
}

/* ------------------------------------------------------------------------ */
/* Incremental insertion                                                    */
/* ------------------------------------------------------------------------ */

struct DetailRow
{
    std::string label;
    size_t      tris;
    double      ms;
    DetailRow(const std::string& l, size_t t, double m) : label(l), tris(t), ms(m) {}
};

static void bench_insert(const Config& c)
{
    section("Incremental insertion",
        "A base triangulation of N uniform random points, followed by " + std::to_string(c.insert_count) +
        " single point insertions. The first insertion builds the internal mesh and is much slower than the "
        "later ones. Break-even is the estimated number of inserted points at which a full rebuild becomes "
        "cheaper: 1 + (build - first insert) / mean(later inserts). Serial = `qdt_execute` + "
        "`qdt_insert_point`; parallel = `qdt_execute_parallel` + `qdt_insert_point_parallel`. Medians of up "
        "to 3 runs.");

    Table t;
    t.col("Base points").lcol("Mode").col("Build (ms)").col("1st insert (ms)").col("Later median (ms)").col(
        "Later mean (ms)").col("Later p95 (ms)").col("Break-even (points)");

    std::vector<DetailRow> detail[2];
    const size_t           m = std::min<size_t>(10, c.insert_count);

    for (size_t k = 0; k < c.insert_sizes.size(); ++k)
    {
        const size_t n = c.insert_sizes[k];
        const bool   last = (k + 1 == c.insert_sizes.size());
        progress("insert, base " + with_commas(n) + " points");

        Rng       rng(c.seed + 31 * n);
        PointList pts = gen_uniform(n, rng);
        PointList extra = gen_uniform(c.insert_count, rng);
        const size_t reps = std::min<size_t>(reps_for(c, n), 3);

        for (int mode = 0; mode < 2; ++mode)
        {
            const bool          par = (mode == 1);
            std::vector<double> builds, firsts, later;

            for (size_t r = 0; r < reps; ++r)
            {
                qdt_t* q = qdt_create();
                const double b0 = now_ms();
                if (par)
                {
                    qdt_execute_parallel(q, pdata(pts), n);
                }
                else
                {
                    qdt_execute(q, pdata(pts), n);
                }
                const double b1 = now_ms();
                const size_t tri0 = qdt_triangles_size(q);
                builds.push_back(b1 - b0);
                if (last && r == 0)
                {
                    detail[mode].push_back(DetailRow(fmt("Most, %zu points left", m), tri0, b1 - b0));
                }

                for (size_t i = 0; i < extra.size(); ++i)
                {
                    const double a = now_ms();
                    const size_t ret = par ? qdt_insert_point_parallel(q, extra[i]) : qdt_insert_point(q, extra[i]);
                    const double b = now_ms();
                    if (ret != n + i + 1)
                    {
                        warn(fmt("insert, base %zu: insertion %zu returned %zu, expected %zu", n, i, ret, n + i + 1));
                    }
                    if (i == 0)
                    {
                        firsts.push_back(b - a);
                    }
                    else
                    {
                        later.push_back(b - a);
                    }
                    if (last && r == 0 && i < 10)
                    {
                        detail[mode].push_back(DetailRow(fmt("insert(%zu)", i), qdt_triangles_size(q), b - a));
                    }
                }
                if (qdt_triangles_size(q) < tri0)
                {
                    warn(fmt("insert, base %zu: triangle count dropped after insertions", n));
                }
                qdt_destroy(q);
            }

            if (last)
            {
                /* "Full": all points, including the first m inserted ones, in one batch call. */
                PointList all = pts;
                all.insert(all.end(), extra.begin(), extra.begin() + static_cast<std::ptrdiff_t>(m));
                size_t       ft = 0;
                const double fm = time_execute(all, par, 1, &ft);
                detail[mode].insert(detail[mode].begin(), DetailRow("Full (" + with_commas(all.size()) + " points)", ft, fm));
            }

            const Stats  sb = summarize(builds);
            const Stats  sf = summarize(firsts);
            const Stats  sl = summarize(later);
            std::string  be = "n/a";
            if (sl.mean > 0.0)
            {
                const double v = 1.0 + (sb.median - sf.median) / sl.mean;
                be = (v < 1.0) ? std::string("<1") : fmt("~%.0f", v);
            }
            const char* mname = par ? "parallel" : "serial";
            t.add(Row()(with_commas(n))(mname)(fmt_ms(sb.median))(fmt_ms(sf.median))(fmt_ms(sl.median))(
                fmt_ms(sl.mean))(fmt_ms(sl.p95))(be));

            rec("insert", mname, n, "build_ms_median", sb.median);
            rec("insert", mname, n, "first_insert_ms_median", sf.median);
            rec("insert", mname, n, "later_insert_ms_median", sl.median);
            rec("insert", mname, n, "later_insert_ms_mean", sl.mean);
            rec("insert", mname, n, "later_insert_ms_p95", sl.p95);
        }
    }
    t.print();

    if (!c.insert_sizes.empty())
    {
        const size_t n = c.insert_sizes.back();
        for (int mode = 0; mode < 2; ++mode)
        {
            subsection(std::string("Per-insert detail, ") + (mode == 1 ? "parallel" : "non-parallel") + ", base " +
                with_commas(n) + " points");
            Table d;
            d.lcol("Operation").col("Triangles").col("Time (ms)");
            for (size_t i = 0; i < detail[mode].size(); ++i)
            {
                d.add(Row()(detail[mode][i].label)(with_commas(detail[mode][i].tris))(fmt_ms(detail[mode][i].ms)));
            }
            d.print();
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Constrained Delaunay triangulation                                       */
/* ------------------------------------------------------------------------ */

/* Bucket grid over the point set, used to pick spatially meaningful constraints. */
struct Grid
{
    int                 g;
    double              cell;
    std::vector<size_t> start; /* g * g + 1 offsets into items */
    std::vector<size_t> items; /* point indices sorted by cell  */

    Grid() : g(0), cell(1.0) {}

    int cidx(double v) const
    {
        int i = static_cast<int>(v / cell);
        if (i < 0) i = 0;
        if (i >= g) i = g - 1;
        return i;
    }
    size_t cell_begin(int cx, int cy) const { return start[static_cast<size_t>(cy) * static_cast<size_t>(g) + static_cast<size_t>(cx)]; }
    size_t cell_end(int cx, int cy) const { return start[static_cast<size_t>(cy) * static_cast<size_t>(g) + static_cast<size_t>(cx) + 1]; }

    void build(const PointList& pts, double extent, double per_cell)
    {
        int gg = static_cast<int>(std::sqrt(static_cast<double>(pts.size()) / per_cell));
        if (gg < 4)
        {
            gg = 4;
        }
        g = gg;
        cell = extent / static_cast<double>(g);
        const size_t cells = static_cast<size_t>(g) * static_cast<size_t>(g);
        start.assign(cells + 1, 0);
        std::vector<size_t> cellid(pts.size());
        for (size_t i = 0; i < pts.size(); ++i)
        {
            cellid[i] = static_cast<size_t>(cidx(pts[i].y)) * static_cast<size_t>(g) + static_cast<size_t>(cidx(pts[i].x));
            ++start[cellid[i] + 1];
        }
        for (size_t k = 0; k < cells; ++k)
        {
            start[k + 1] += start[k];
        }
        items.assign(pts.size(), 0);
        std::vector<size_t> fill(start.begin(), start.end() - 1);
        for (size_t i = 0; i < pts.size(); ++i)
        {
            items[fill[cellid[i]]++] = i;
        }
    }
};

static qdt_edge make_edge(size_t a, size_t b)
{
    qdt_edge e;
    e.v0 = a;
    e.v1 = b;
    return e;
}

/*
 * Short constraints: one segment per 3x3 block of grid cells, from a point of
 * the centre cell to another point of the surrounding 3x3 block. Blocks are
 * three cells apart, so no two segments can cross or touch.
 */
static std::vector<qdt_edge> make_short_constraints(const Grid& gr, Rng& rng, size_t want)
{
    std::vector<qdt_edge> out;
    std::vector<size_t>   cand;
    for (int cy = 0; cy < gr.g; cy += 3)
    {
        for (int cx = 0; cx < gr.g; cx += 3)
        {
            const size_t b0 = gr.cell_begin(cx, cy), b1 = gr.cell_end(cx, cy);
            if (b0 == b1)
            {
                continue;
            }
            const size_t a = gr.items[b0 + rng.index(b1 - b0)];
            cand.clear();
            for (int dy = -1; dy <= 1; ++dy)
            {
                for (int dx = -1; dx <= 1; ++dx)
                {
                    const int x = cx + dx, y = cy + dy;
                    if (x < 0 || y < 0 || x >= gr.g || y >= gr.g)
                    {
                        continue;
                    }
                    for (size_t k = gr.cell_begin(x, y); k < gr.cell_end(x, y); ++k)
                    {
                        if (gr.items[k] != a)
                        {
                            cand.push_back(gr.items[k]);
                        }
                    }
                }
            }
            if (!cand.empty())
            {
                out.push_back(make_edge(a, cand[rng.index(cand.size())]));
            }
        }
    }
    shuffle_vec(out, rng);
    if (out.size() > want)
    {
        out.resize(want);
    }
    return out;
}

/*
 * Long chords: in evenly spaced grid rows, one segment from the leftmost to
 * the rightmost point of the row. Both end points lie in the same row band, so
 * chords of different rows never cross, and each one spans the whole domain
 * (and therefore crosses about sqrt(n) triangles).
 */
static std::vector<qdt_edge> make_long_chords(const Grid& gr, size_t want)
{
    std::vector<qdt_edge> out;
    size_t                L = want;
    if (L > static_cast<size_t>(gr.g) / 2)
    {
        L = static_cast<size_t>(gr.g) / 2;
    }
    for (size_t k = 0; k < L; ++k)
    {
        const int row = static_cast<int>((static_cast<double>(k) + 0.5) * static_cast<double>(gr.g) / static_cast<double>(L));
        int       lx = -1, rx = -1;
        for (int x = 0; x < gr.g; ++x)
        {
            if (gr.cell_begin(x, row) != gr.cell_end(x, row))
            {
                lx = x;
                break;
            }
        }
        for (int x = gr.g - 1; x >= 0; --x)
        {
            if (gr.cell_begin(x, row) != gr.cell_end(x, row))
            {
                rx = x;
                break;
            }
        }
        if (lx < 0 || rx <= lx)
        {
            continue;
        }
        out.push_back(make_edge(gr.items[gr.cell_begin(lx, row)], gr.items[gr.cell_begin(rx, row)]));
    }
    return out;
}

/*
 * Closed polygon: all points inside a thin annulus around the centre of the
 * domain, sorted by angle and joined consecutively. The result is a simple
 * (star-shaped) polygon, so the ring edges never cross.
 */
static std::vector<qdt_edge> make_ring(const PointList& pts, size_t target)
{
    std::vector<qdt_edge> out;
    const double          cx = kExtent * 0.5, cy = kExtent * 0.5, r0 = kExtent * 0.3;
    const double          density = static_cast<double>(pts.size()) / (kExtent * kExtent);
    double                width = static_cast<double>(target) / (density * 2.0 * kPi * r0);
    if (width > kExtent * 0.1)
    {
        width = kExtent * 0.1;
    }

    std::vector<std::pair<double, size_t>> v;
    for (size_t i = 0; i < pts.size(); ++i)
    {
        const double dx = pts[i].x - cx, dy = pts[i].y - cy;
        const double d = std::sqrt(dx * dx + dy * dy);
        if (std::fabs(d - r0) <= width * 0.5)
        {
            v.push_back(std::make_pair(std::atan2(dy, dx), i));
        }
    }
    if (v.size() < 8)
    {
        return out;
    }
    std::sort(v.begin(), v.end());
    for (size_t k = 0; k < v.size(); ++k)
    {
        out.push_back(make_edge(v[k].second, v[(k + 1) % v.size()].second));
    }
    return out;
}

struct ConsResult
{
    double first_ms; /* first qdt_insert_constraint call: includes the one-time mesh build */
    double rest_ms;  /* all other constraints */
    size_t accepted;
    size_t tris_before;
    size_t tris_after;
};

/* Resets the context with qdt_execute (untimed), then inserts the constraints. */
static ConsResult run_constraints(qdt_t* q, const PointList& pts, const std::vector<qdt_edge>& es, bool batch)
{
    qdt_execute(q, pdata(pts), pts.size());
    ConsResult r;
    r.tris_before = qdt_triangles_size(q);

    const double t0 = now_ms();
    size_t       acc = qdt_insert_constraint(q, es[0].v0, es[0].v1);
    const double t1 = now_ms();
    if (batch)
    {
        acc += qdt_insert_constraints_batch(q, &es[1], es.size() - 1);
    }
    else
    {
        for (size_t i = 1; i < es.size(); ++i)
        {
            acc += qdt_insert_constraint(q, es[i].v0, es[i].v1);
        }
    }
    const double t2 = now_ms();

    r.first_ms = t1 - t0;
    r.rest_ms = t2 - t1;
    r.accepted = acc;
    r.tris_after = qdt_triangles_size(q);
    return r;
}

struct ConsSet
{
    std::string           name;
    std::vector<qdt_edge> edges;
};

static void bench_constraint(const Config& c)
{
    section("Constrained Delaunay triangulation",
        "A base triangulation of N uniform random points, then constraints from three sets: \"short\" "
        "(segments between nearby points, none of them crossing; some are already Delaunay edges), \"long "
        "chords\" (segments spanning the whole domain, each crossing about sqrt(N) triangles) and \"ring\" (a "
        "closed polygon through the points of a thin annulus). The first constraint is inserted with "
        "`qdt_insert_constraint` and includes the one-time internal mesh build; the others are inserted either "
        "with one `qdt_insert_constraints_batch` call or with a loop of `qdt_insert_constraint` calls. Every "
        "constraint of these sets is valid and must be accepted. Medians of up to 3 runs.");

    Table t;
    t.col("Points").lcol("Constraint set").col("Constraints").col("1st call + mesh build (ms)").col("Batch (ms)").col(
        "Batch (us/constraint)").col("Single calls (ms)").col("Single (us/constraint)").col("Accepted");

    for (size_t k = 0; k < c.constraint_sizes.size(); ++k)
    {
        const size_t n = c.constraint_sizes[k];
        progress("constraints, " + with_commas(n) + " points");
        Rng       rng(c.seed + 77 * n);
        PointList pts = gen_uniform(n, rng);
        Grid      gr;
        gr.build(pts, kExtent, 2.0);

        std::vector<ConsSet> sets(3);
        sets[0].name = "short";
        sets[0].edges = make_short_constraints(gr, rng, c.short_constraints);
        sets[1].name = "long chords";
        sets[1].edges = make_long_chords(gr, c.long_chords);
        sets[2].name = "ring";
        sets[2].edges = make_ring(pts, c.ring_points);

        const size_t reps = std::min<size_t>(reps_for(c, n), 3);
        qdt_t* q = qdt_create();

        for (size_t s = 0; s < sets.size(); ++s)
        {
            const std::vector<qdt_edge>& es = sets[s].edges;
            if (es.size() < 2)
            {
                warn(fmt("constraints, %zu points, set '%s': too few constraints could be generated (%zu)", n,
                    sets[s].name.c_str(), es.size()));
                continue;
            }

            std::vector<double> first, brest, srest;
            size_t              acc_b = 0, acc_s = 0;
            for (size_t r = 0; r < reps; ++r)
            {
                const ConsResult b = run_constraints(q, pts, es, true);
                const ConsResult x = run_constraints(q, pts, es, false);
                first.push_back(b.first_ms);
                brest.push_back(b.rest_ms);
                srest.push_back(x.rest_ms);
                acc_b = b.accepted;
                acc_s = x.accepted;
                if (b.tris_before != b.tris_after || x.tris_before != x.tris_after)
                {
                    warn(fmt("constraints, %zu points, set '%s': triangle count changed (%zu -> %zu)", n,
                        sets[s].name.c_str(), b.tris_before, b.tris_after));
                }
            }
            if (acc_b != es.size() || acc_s != es.size())
            {
                warn(fmt("constraints, %zu points, set '%s': accepted %zu (batch) / %zu (single) of %zu", n,
                    sets[s].name.c_str(), acc_b, acc_s, es.size()));
            }

            const double f = summarize(first).median;
            const double bm = summarize(brest).median;
            const double sm = summarize(srest).median;
            const double nr = static_cast<double>(es.size() - 1);
            t.add(Row()(with_commas(n))(sets[s].name)(with_commas(es.size()))(fmt_ms(f))(fmt_ms(bm))(
                fmt_us(bm * 1000.0 / nr))(fmt_ms(sm))(fmt_us(sm * 1000.0 / nr))(
                    with_commas(acc_b) + " / " + with_commas(acc_s)));

            rec("constraint", sets[s].name, n, "constraints", static_cast<double>(es.size()));
            rec("constraint", sets[s].name, n, "first_call_ms_median", f);
            rec("constraint", sets[s].name, n, "batch_ms_median", bm);
            rec("constraint", sets[s].name, n, "single_ms_median", sm);
        }
        qdt_destroy(q);
    }
    t.print();
    std::printf("Accepted is shown as batch / single-call counts.\n\n");
}

/* ------------------------------------------------------------------------ */
/* Reading results                                                          */
/* ------------------------------------------------------------------------ */

static void bench_query(const Config& c)
{
    const size_t n = c.query_n;
    section("Reading the result", "`qdt_get_triangle` over all triangles of a uniform random point set.");
    if (n == 0)
    {
        std::printf("(skipped: --max-points is too small)\n\n");
        return;
    }
    progress("query, " + with_commas(n) + " points");

    Rng       rng(c.seed + 555);
    PointList pts = gen_uniform(n, rng);
    qdt_t* q = qdt_create();
    qdt_execute(q, pdata(pts), n);
    const size_t tris = qdt_triangles_size(q);

    std::vector<double> t;
    const size_t        reps = reps_for(c, n);
    for (size_t r = 0; r < reps + 1; ++r) /* the first pass is a warm-up */
    {
        uint64_t     sum = 0;
        const double t0 = now_ms();
        for (size_t i = 0; i < tris; ++i)
        {
            size_t T[3];
            if (qdt_get_triangle(q, T, i))
            {
                sum += T[0] ^ T[1] ^ T[2];
            }
        }
        const double t1 = now_ms();
        g_sink = g_sink + sum;
        if (r > 0)
        {
            t.push_back(t1 - t0);
        }
    }
    qdt_destroy(q);

    const double med = summarize(t).median;
    Table        tb;
    tb.col("Points").col("Triangles").col("Time (ms)").col("ns per call").col("Mcalls/s");
    tb.add(Row()(with_commas(n))(with_commas(tris))(fmt_ms(med))(fmt("%.1f", med * 1e6 / static_cast<double>(tris)))(
        fmt("%.1f", static_cast<double>(tris) / med / 1000.0)));
    tb.print();
    rec("query", "get_triangle", n, "ms_median", med);
}

/* ------------------------------------------------------------------------ */
/* Independent contexts on several threads                                  */
/* ------------------------------------------------------------------------ */

static double run_concurrent_once(const PointList& pts, unsigned threads, size_t* tris_out, bool* consistent)
{
    std::vector<size_t>      tris(threads, 0);
    std::vector<std::thread> th;
    th.reserve(threads);

    const double t0 = now_ms();
    for (unsigned i = 0; i < threads; ++i)
    {
        th.push_back(std::thread([&pts, &tris, i]() {
            qdt_t* q = qdt_create();
            qdt_execute(q, pdata(pts), pts.size());
            tris[i] = qdt_triangles_size(q);
            qdt_destroy(q);
            }));
    }
    for (unsigned i = 0; i < threads; ++i)
    {
        th[i].join();
    }
    const double t1 = now_ms();

    *consistent = true;
    for (unsigned i = 1; i < threads; ++i)
    {
        if (tris[i] != tris[0])
        {
            *consistent = false;
        }
    }
    *tris_out = tris[0];
    return t1 - t0;
}

static void bench_concurrent(const Config& c)
{
    const size_t n = c.concurrent_n;
    unsigned     hw = std::thread::hardware_concurrency();
    if (hw == 0)
    {
        hw = 1;
    }
    section("Independent contexts on several threads",
        "Every thread triangulates the same " + with_commas(n) +
        " uniform random points with `qdt_execute` on its own context (the contexts share no state). "
        "Scaling is the aggregate throughput divided by the throughput of one thread; efficiency is "
        "scaling divided by the thread count. Hardware threads: " +
        std::to_string(hw) + ". Medians of " + std::to_string(c.reps) + " runs.");
    if (n == 0)
    {
        std::printf("(skipped: --max-points is too small)\n\n");
        return;
    }

    Rng       rng(c.seed + 666);
    PointList pts = gen_uniform(n, rng);

    std::vector<unsigned> counts;
    for (unsigned th = 1; th < hw; th *= 2)
    {
        counts.push_back(th);
    }
    counts.push_back(hw);

    Table  t;
    double base_throughput = 0.0;
    size_t ref_tris = 0;
    t.col("Threads").col("Wall time (ms)").col("Aggregate (Mpts/s)").col("Scaling").col("Efficiency");
    for (size_t k = 0; k < counts.size(); ++k)
    {
        const unsigned th = counts[k];
        progress("concurrent, " + std::to_string(th) + " threads");

        std::vector<double> times;
        for (size_t r = 0; r < c.reps + 1; ++r) /* the first run is a warm-up */
        {
            size_t tris = 0;
            bool   same = true;
            const double ms = run_concurrent_once(pts, th, &tris, &same);
            if (!same)
            {
                warn(fmt("concurrent, %u threads: contexts produced different triangle counts", th));
            }
            if (ref_tris == 0)
            {
                ref_tris = tris;
            }
            else if (tris != ref_tris)
            {
                warn(fmt("concurrent, %u threads: triangle count %zu differs from %zu", th, tris, ref_tris));
            }
            if (r > 0)
            {
                times.push_back(ms);
            }
        }
        const double med = summarize(times).median;
        const double throughput = static_cast<double>(th) * static_cast<double>(n) / med / 1000.0; /* Mpts/s */
        if (k == 0)
        {
            base_throughput = throughput;
        }
        const double scaling = throughput / base_throughput;
        t.add(Row()(std::to_string(th))(fmt_ms(med))(fmt("%.2f", throughput))(fmt_x(scaling))(
            fmt("%.0f%%", 100.0 * scaling / static_cast<double>(th))));
        rec("concurrent", std::to_string(th) + " threads", n, "wall_ms_median", med);
        rec("concurrent", std::to_string(th) + " threads", n, "mpts_per_s", throughput);
    }
    t.print();
}

/* ------------------------------------------------------------------------ */
/* System information                                                       */
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
    if (reg_string(HKEY_LOCAL_MACHINE, "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0", "ProcessorNameString", v) &&
        !v.empty())
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

/* Number of physical cores, 0 if unknown. */
static unsigned physical_cores()
{
#if defined(_WIN32)
    DWORD len = 0;
    GetLogicalProcessorInformation(NULL, &len);
    if (GetLastError() != ERROR_INSUFFICIENT_BUFFER || len == 0)
    {
        return 0;
    }
    std::vector<SYSTEM_LOGICAL_PROCESSOR_INFORMATION> buf(len / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION));
    if (buf.empty() || !GetLogicalProcessorInformation(&buf[0], &len))
    {
        return 0;
    }
    unsigned cores = 0;
    for (size_t i = 0; i < buf.size(); ++i)
    {
        if (buf[i].Relationship == RelationProcessorCore)
        {
            ++cores;
        }
    }
    return cores;
#else
    FILE* f = std::fopen("/proc/cpuinfo", "r");
    if (!f)
    {
        return 0;
    }
    std::set<std::pair<long, long>> cores;
    long                            phys = -1, core = -1;
    char                            line[512];
    while (std::fgets(line, sizeof(line), f))
    {
        if (line[0] == '\n' || line[0] == '\r')
        {
            if (core >= 0)
            {
                cores.insert(std::make_pair(phys, core));
            }
            phys = -1;
            core = -1;
            continue;
        }
        const char* colon = std::strchr(line, ':');
        if (!colon)
        {
            continue;
        }
        if (std::strncmp(line, "physical id", 11) == 0)
        {
            phys = std::atol(colon + 1);
        }
        else if (std::strncmp(line, "core id", 7) == 0)
        {
            core = std::atol(colon + 1);
        }
    }
    if (core >= 0)
    {
        cores.insert(std::make_pair(phys, core));
    }
    std::fclose(f);
    return static_cast<unsigned>(cores.size());
#endif
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
    bool         have = false;
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
            name = fmt("Windows %lu.%lu", static_cast<unsigned long>(vi.dwMajorVersion),
                static_cast<unsigned long>(vi.dwMinorVersion));
        }
    }
    std::string display;
    if (reg_string(HKEY_LOCAL_MACHINE, "SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion", "DisplayVersion", display) &&
        !display.empty())
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
#elif defined(QDT_BENCH_LINUX)
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

/* ------------------------------------------------------------------------ */
/* Memory use                                                               */
/* ------------------------------------------------------------------------ */

/*
 * Every measurement runs in a process of its own so that it starts from a
 * clean heap and the figures do not depend on the scenarios that ran before.
 * mem_measure() fills out[5]:
 *   out[0]: bit 0 = values valid, bit 1 = batch peak valid, bit 2 = mesh peak valid
 *   out[1]: bytes retained after the batch call
 *   out[2]: peak bytes above the baseline during the batch call
 *   out[3]: bytes added by the first insertion (internal mesh)
 *   out[4]: peak bytes above the post-batch level during that insertion
 */

static uint64_t bytes_delta(uint64_t after, uint64_t before)
{
    return after > before ? after - before : 0ULL;
}

#if defined(QDT_BENCH_LINUX)

static bool read_kb(const char* key, uint64_t* out)
{
    std::string s;
    if (!proc_value("/proc/self/status", key, s))
    {
        return false;
    }
    *out = static_cast<uint64_t>(std::strtoull(s.c_str(), NULL, 10));
    return true;
}

/* Resets VmHWM (peak RSS) to the current RSS. Needs Linux >= 4.0 and write access to /proc/self/clear_refs. */
static bool reset_peak_rss()
{
    FILE* f = std::fopen("/proc/self/clear_refs", "w");
    if (!f)
    {
        return false;
    }
    const bool wrote = std::fputs("5", f) >= 0;
    const bool closed = std::fclose(f) == 0;
    return wrote && closed;
}

/*
 * Gives memory that was freed but is still resident back to the operating
 * system. Without it a forked child inherits the heap of the parent, which is
 * full of freed-but-resident pages left by earlier scenarios; the library
 * would then reuse those pages and the RSS growth would be badly
 * underestimated. Only glibc provides malloc_trim.
 */
static void trim_heap()
{
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
}

/*
 * Runs in a forked child. The heap is trimmed before the baseline and before
 * every RSS reading, so "retained" approximates the live memory of the
 * context; the peak values are read before trimming.
 */
static void mem_child(const PointList& pts, bool parallel, int wfd)
{
    uint64_t out[5] = { 0, 0, 0, 0, 0 };
    uint64_t rss0 = 0, rss1 = 0, rss2 = 0, hwm1 = 0, hwm2 = 0;

    trim_heap();
    bool peak = reset_peak_rss();
    bool ok = read_kb("VmRSS", &rss0);

    qdt_t* q = qdt_create();
    if (parallel)
    {
        qdt_execute_parallel(q, pdata(pts), pts.size());
    }
    else
    {
        qdt_execute(q, pdata(pts), pts.size());
    }
    ok = ok && read_kb("VmHWM", &hwm1);
    trim_heap();
    ok = ok && read_kb("VmRSS", &rss1);

    if (peak)
    {
        peak = reset_peak_rss();
    }
    const qdt_vec2d extra = P(kExtent * 0.5123, kExtent * 0.4877);
    if (parallel)
    {
        qdt_insert_point_parallel(q, extra);
    }
    else
    {
        qdt_insert_point(q, extra);
    }
    ok = ok && read_kb("VmHWM", &hwm2);
    trim_heap();
    ok = ok && read_kb("VmRSS", &rss2);

    out[0] = ok ? (peak ? 7 : 1) : 0;
    out[1] = bytes_delta(rss1, rss0) * 1024ULL;
    out[2] = peak ? bytes_delta(hwm1, rss0) * 1024ULL : 0;
    out[3] = bytes_delta(rss2, rss1) * 1024ULL;
    out[4] = peak ? bytes_delta(hwm2, rss1) * 1024ULL : 0;

    const ssize_t w = write(wfd, out, sizeof(out));
    (void)w;
}

static bool mem_measure(size_t n, uint64_t seed, bool parallel, uint64_t out[5])
{
    Rng       rng(seed);
    PointList pts = gen_uniform(n, rng);

    int fd[2];
    if (pipe(fd) != 0)
    {
        return false;
    }
    std::fflush(stdout);
    std::fflush(stderr);

    const pid_t pid = fork();
    if (pid < 0)
    {
        close(fd[0]);
        close(fd[1]);
        return false;
    }
    if (pid == 0)
    {
        close(fd[0]);
        try
        {
            mem_child(pts, parallel, fd[1]);
        }
        catch (...)
        {
            _exit(1);
        }
        close(fd[1]);
        _exit(0);
    }

    close(fd[1]);
    const size_t want = 5 * sizeof(uint64_t);
    size_t       got = 0;
    char* buf = reinterpret_cast<char*>(out);
    while (got < want)
    {
        const ssize_t r = read(fd[0], buf + got, want - got);
        if (r <= 0)
        {
            break;
        }
        got += static_cast<size_t>(r);
    }
    close(fd[0]);
    int status = 0;
    waitpid(pid, &status, 0);
    return got == want && (out[0] & 1) != 0;
}

static const char* const kMemNote =
"Resident set size (VmRSS) measured around the library calls, in a freshly forked process per measurement; "
"the input point array is part of the baseline. Freed memory is returned to the operating system "
"(`malloc_trim`, glibc only) before each reading, so \"retained\" approximates the live memory of the context.";

#elif defined(_WIN32)

static bool win_commit(uint64_t* commit, uint64_t* peak_commit)
{
    PROCESS_MEMORY_COUNTERS_EX pmc;
    std::memset(&pmc, 0, sizeof(pmc));
    pmc.cb = sizeof(pmc);
    if (!GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc),
        static_cast<DWORD>(sizeof(pmc))))
    {
        return false;
    }
    *commit = static_cast<uint64_t>(pmc.PrivateUsage);
    *peak_commit = static_cast<uint64_t>(pmc.PeakPagefileUsage);
    return true;
}

/*
 * Entry point of the child process (`qdt_bench --mem-child N SEED PARALLEL`).
 * Windows has no fork and no way to reset the peak counters, so the parent
 * starts a fresh copy of this program for every measurement. The child
 * generates the points, takes the baseline, runs the batch call and the first
 * insertion, and prints one result line. The peak commit of a fresh process
 * is only attributable to a call if the call raised it above the value it had
 * before the call; otherwise the peak is reported as unavailable.
 */
static int mem_child_main(int argc, char** argv)
{
    if (argc < 5)
    {
        return 2;
    }
    const size_t   n = static_cast<size_t>(std::strtoull(argv[2], NULL, 10));
    const uint64_t seed = static_cast<uint64_t>(std::strtoull(argv[3], NULL, 10));
    const bool     par = std::atoi(argv[4]) != 0;

    Rng       rng(seed);
    PointList pts = gen_uniform(n, rng);

    uint64_t c0 = 0, p0 = 0, c1 = 0, p1 = 0, c2 = 0, p2 = 0;
    bool     ok = win_commit(&c0, &p0);

    qdt_t* q = qdt_create();
    if (par)
    {
        qdt_execute_parallel(q, pdata(pts), pts.size());
    }
    else
    {
        qdt_execute(q, pdata(pts), pts.size());
    }
    ok = ok && win_commit(&c1, &p1);

    const qdt_vec2d extra = P(kExtent * 0.5123, kExtent * 0.4877);
    if (par)
    {
        qdt_insert_point_parallel(q, extra);
    }
    else
    {
        qdt_insert_point(q, extra);
    }
    ok = ok && win_commit(&c2, &p2);
    qdt_destroy(q);

    /*
     * The result of a triangulation always occupies memory, so the commit
     * charge must have grown after the batch call and after the first
     * insertion. If it did not, the system does not provide meaningful
     * counters (Wine, for example, reports PrivateUsage as 0): report a
     * failure instead of printing zeros.
     */
    if (ok && (c1 <= c0 || c2 <= c1))
    {
        std::fprintf(stderr, "  [bench] the process memory counters do not change (PrivateUsage unavailable?)\n");
        ok = false;
    }

    uint64_t flags = ok ? 1 : 0;
    if (ok && p1 > p0)
    {
        flags |= 2;
    }
    if (ok && p2 > p1)
    {
        flags |= 4;
    }
    std::printf("QDTMEM %llu %llu %llu %llu %llu\n", static_cast<unsigned long long>(flags),
        static_cast<unsigned long long>(bytes_delta(c1, c0)),
        static_cast<unsigned long long>((flags & 2) ? bytes_delta(p1, c0) : 0),
        static_cast<unsigned long long>(bytes_delta(c2, c1)),
        static_cast<unsigned long long>((flags & 4) ? bytes_delta(p2, c1) : 0));
    std::fflush(stdout);
    return ok ? 0 : 1;
}

/* Runs a command line, collects its standard output, returns true if it exited with code 0. */
static bool win_run_capture(const std::string& cmdline, std::string& output)
{
    SECURITY_ATTRIBUTES sa;
    std::memset(&sa, 0, sizeof(sa));
    sa.nLength = sizeof(sa);
    sa.bInheritHandle = TRUE;

    HANDLE rd = NULL, wr = NULL;
    if (!CreatePipe(&rd, &wr, &sa, 0))
    {
        return false;
    }
    SetHandleInformation(rd, HANDLE_FLAG_INHERIT, 0);

    STARTUPINFOA si;
    std::memset(&si, 0, sizeof(si));
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdOutput = wr;
    si.hStdError = GetStdHandle(STD_ERROR_HANDLE);
    si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);

    PROCESS_INFORMATION pi;
    std::memset(&pi, 0, sizeof(pi));
    std::vector<char> cl(cmdline.begin(), cmdline.end());
    cl.push_back('\0'); /* CreateProcessA may modify the buffer */

    const BOOL started = CreateProcessA(NULL, &cl[0], NULL, NULL, TRUE, 0, NULL, NULL, &si, &pi);
    CloseHandle(wr);
    if (!started)
    {
        CloseHandle(rd);
        return false;
    }

    char  buf[256];
    DWORD got = 0;
    while (ReadFile(rd, buf, static_cast<DWORD>(sizeof(buf)), &got, NULL) && got > 0)
    {
        output.append(buf, static_cast<size_t>(got));
    }
    CloseHandle(rd);

    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD code = 1;
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    return code == 0;
}

static bool mem_measure(size_t n, uint64_t seed, bool parallel, uint64_t out[5])
{
    std::vector<char> path(32768, 0);
    const DWORD       len = GetModuleFileNameA(NULL, &path[0], static_cast<DWORD>(path.size()));
    if (len == 0 || len >= path.size())
    {
        return false;
    }
    const std::string cmd = "\"" + std::string(&path[0], static_cast<size_t>(len)) + "\" --mem-child " +
        std::to_string(static_cast<unsigned long long>(n)) + " " +
        std::to_string(static_cast<unsigned long long>(seed)) + " " + (parallel ? "1" : "0");

    std::string text;
    if (!win_run_capture(cmd, text))
    {
        return false;
    }
    const size_t pos = text.find("QDTMEM");
    if (pos == std::string::npos)
    {
        return false;
    }
    unsigned long long v[5] = { 0, 0, 0, 0, 0 };
    if (std::sscanf(text.c_str() + pos, "QDTMEM %llu %llu %llu %llu %llu", &v[0], &v[1], &v[2], &v[3], &v[4]) != 5)
    {
        return false;
    }
    for (int i = 0; i < 5; ++i)
    {
        out[i] = static_cast<uint64_t>(v[i]);
    }
    return (out[0] & 1) != 0;
}

static const char* const kMemNote =
"Committed private bytes (PrivateUsage) measured around the library calls, in a freshly started process per "
"measurement; the input point array is part of the baseline.";

#endif

static void bench_memory(const Config& c)
{
#if defined(QDT_BENCH_HAVE_MEMORY)
    section("Memory use",
        std::string(kMemNote) +
        " \"Retained\" is the growth after the call. \"Peak\" is the highest level during the call minus the "
        "level before it. Mesh = the internal half-edge mesh built by the first insertion after a batch "
        "call. Bytes are per input point.");

    Table t;
    bool  any_na = false;
    t.col("Points").lcol("Mode").col("Retained after batch (B/pt)").col("Peak during batch (B/pt)").col(
        "Mesh retained (B/pt)").col("Peak during mesh build (B/pt)").col("Retained, batch + mesh (MiB)");

    for (size_t k = 0; k < c.memory_sizes.size(); ++k)
    {
        const size_t n = c.memory_sizes[k];
        progress("memory, " + with_commas(n) + " points");

        for (int mode = 0; mode < 2; ++mode)
        {
            const bool par = (mode == 1);
            uint64_t   v[5] = { 0, 0, 0, 0, 0 };
            if (!mem_measure(n, c.seed + 99 * n, par, v))
            {
                warn(fmt("memory, %zu points, %s: measurement failed (the process memory counters may be "
                    "unavailable, e.g. under Wine)",
                    n, par ? "parallel" : "serial"));
                continue;
            }
            const double dn = static_cast<double>(n);
            const bool   peak_batch = (v[0] & 2) != 0;
            const bool   peak_mesh = (v[0] & 4) != 0;
            if (!peak_batch || !peak_mesh)
            {
                any_na = true;
            }
            t.add(Row()(with_commas(n))(par ? "parallel" : "serial")(fmt("%.1f", static_cast<double>(v[1]) / dn))(
                peak_batch ? fmt("%.1f", static_cast<double>(v[2]) / dn) : std::string("n/a"))(
                    fmt("%.1f", static_cast<double>(v[3]) / dn))(
                        peak_mesh ? fmt("%.1f", static_cast<double>(v[4]) / dn) : std::string("n/a"))(
                            fmt("%.1f", static_cast<double>(v[1] + v[3]) / (1024.0 * 1024.0))));

            const char* mname = par ? "parallel" : "serial";
            rec("memory", mname, n, "retained_batch_bytes", static_cast<double>(v[1]));
            rec("memory", mname, n, "peak_batch_bytes", static_cast<double>(v[2]));
            rec("memory", mname, n, "retained_mesh_bytes", static_cast<double>(v[3]));
            rec("memory", mname, n, "peak_mesh_bytes", static_cast<double>(v[4]));
        }
    }
    t.print();
    if (any_na)
    {
        std::printf("n/a: the peak could not be determined on this system for that call.\n\n");
    }
#else
    (void)c;
    section("Memory use", "Skipped: the memory scenario needs Linux or Windows.");
#endif
}

/* ------------------------------------------------------------------------ */
/* Optional comparison with CGAL                                            */
/* ------------------------------------------------------------------------ */

#if defined(QDT_BENCH_CGAL)

typedef CGAL::Exact_predicates_inexact_constructions_kernel CgalK;
typedef CGAL::Delaunay_triangulation_2<CgalK>               CgalDT;

static void bench_cgal(const Config& c)
{
    std::string ver = "unknown version";
#if defined(CGAL_VERSION_STR)
    ver = CGAL_VERSION_STR;
#endif
    section("Comparison with CGAL",
        "CGAL " + ver +
        " `Delaunay_triangulation_2` with `Exact_predicates_inexact_constructions_kernel`, inserting the "
        "whole point range at once (CGAL sorts spatially itself). Point conversion is not timed. Uniform "
        "random points, medians of " + std::to_string(c.reps) + " runs (3 for the largest sizes). The speedup "
        "columns are CGAL time divided by QDT time. Make sure CGAL was compiled with NDEBUG.");

    Table t;
    t.col("Points").col("CGAL (ms)").col("QDT serial (ms)").col("QDT parallel (ms)").col("Serial speedup").col(
        "Parallel speedup").col("Triangles (CGAL)").col("Triangles (QDT)");

    for (size_t k = 0; k < c.cgal_sizes.size(); ++k)
    {
        const size_t n = c.cgal_sizes[k];
        progress("cgal, " + with_commas(n) + " points");
        Rng       rng(c.seed + n);
        PointList pts = gen_uniform(n, rng);

        std::vector<CgalK::Point_2> cp;
        cp.reserve(n);
        for (size_t i = 0; i < n; ++i)
        {
            cp.push_back(CgalK::Point_2(pts[i].x, pts[i].y));
        }

        const size_t reps = reps_for(c, n);
        const bool   warm = (n <= 1000000);
        size_t       faces = 0;
        std::vector<double> ct;
        for (size_t r = 0; r < reps + (warm ? 1 : 0); ++r)
        {
            CgalDT       dt;
            const double t0 = now_ms();
            dt.insert(cp.begin(), cp.end());
            const double t1 = now_ms();
            faces = static_cast<size_t>(dt.number_of_faces());
            if (!(warm && r == 0))
            {
                ct.push_back(t1 - t0);
            }
        }

        const BatchResult s = measure_batch(pts, false, reps);
        const BatchResult p = measure_batch(pts, true, reps);
        if (s.tris != faces)
        {
            warn(fmt("cgal, %zu points: CGAL has %zu triangles, QDT has %zu", n, faces, s.tris));
        }

        const double cm = summarize(ct).median;
        t.add(Row()(with_commas(n))(fmt_ms(cm))(fmt_ms(s.ms.median))(fmt_ms(p.ms.median))(fmt_x(cm / s.ms.median))(
            fmt_x(cm / p.ms.median))(with_commas(faces))(with_commas(s.tris)));

        rec("cgal", "cgal", n, "ms_median", cm);
        rec("cgal", "qdt_serial", n, "ms_median", s.ms.median);
        rec("cgal", "qdt_parallel", n, "ms_median", p.ms.median);
    }
    t.print();
}

#endif

/* ------------------------------------------------------------------------ */
/* Environment                                                              */
/* ------------------------------------------------------------------------ */

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

static bool built_optimized()
{
#if defined(__GNUC__) || defined(__clang__)
#if defined(__OPTIMIZE__)
    return true;
#else
    return false;
#endif
#elif defined(_MSC_VER)
#if defined(_DEBUG)
    return false;
#else
    return true;
#endif
#else
    return true;
#endif
}

static void print_environment(const Config& c)
{
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0)
    {
        hw = 1;
    }
    const unsigned cores = physical_cores();
    const unsigned long long mem = total_memory_bytes();

    char        when[64];
    std::time_t now = std::time(NULL);
    std::strftime(when, sizeof(when), "%Y-%m-%d %H:%M UTC", std::gmtime(&now));

    std::printf("## Performance\n\n");
    std::printf("The figures below were produced by `test/qdt_bench.cpp`. Timings are wall-clock milliseconds "
        "around the library calls only (point generation and the creation and destruction of contexts "
        "are not included); tables show medians. Absolute numbers depend on the CPU, compiler, compiler "
        "flags and thread count, so the relative comparisons are more useful than the absolute values.\n\n");
    std::printf("### Test environment\n\n");
    std::printf("- Date: %s\n", when);
    std::printf("- CPU: %s\n", cpu_name().c_str());
    if (cores != 0)
    {
        std::printf("- Cores / hardware threads: %u / %u (`qdt_execute_parallel` uses min(threads, 64))\n", cores, hw);
    }
    else
    {
        std::printf("- Hardware threads: %u (`qdt_execute_parallel` uses min(threads, 64))\n", hw);
    }
    if (mem != 0)
    {
        std::printf("- Memory: %.1f GiB\n", static_cast<double>(mem) / (1024.0 * 1024.0 * 1024.0));
    }
    else
    {
        std::printf("- Memory: unknown\n");
    }
    std::printf("- OS: %s\n", os_name().c_str());
    std::printf("- Compiler (benchmark program): %s, %zu-bit\n", compiler_string().c_str(), sizeof(void*) * 8);
    std::printf("- Parameters: repeat = %zu, seed = %llu%s\n", c.reps, static_cast<unsigned long long>(c.seed),
        c.quick ? ", quick mode" : "");
    if (c.max_points != 0)
    {
        std::printf("- Sizes limited to at most %s points\n", with_commas(c.max_points).c_str());
    }
    for (size_t i = 0; i < c.notes.size(); ++i)
    {
        std::printf("- %s\n", c.notes[i].c_str());
    }
    if (!built_optimized())
    {
        std::printf("\n> **Warning:** this benchmark was compiled without optimization. Only the benchmark's own "
            "flags can be checked; libqdt must be built with optimizations as well. Numbers from "
            "unoptimized builds are meaningless.\n");
        warn("the benchmark program was compiled without optimization");
    }
    std::printf("\n");
    std::fflush(stdout);
}

/* ------------------------------------------------------------------------ */
/* Driver                                                                   */
/* ------------------------------------------------------------------------ */

static const char* const kScenarioNames[] = { "batch", "dist", "insert", "constraint", "query", "concurrent", "memory", "cgal" };

static void usage(const char* argv0)
{
    std::fprintf(stderr,
        "usage: %s [options]\n"
        "  --output FILE      write the Markdown report to FILE (UTF-8, LF) instead of stdout\n"
        "  --note TEXT        add a line to the \"Test environment\" list (repeatable)\n"
        "  --quick            small sizes only (smoke run)\n"
        "  --repeat N         repetitions per measurement (default 5)\n"
        "  --seed N           seed of the pseudo random generator\n"
        "  --scenario LIST    comma separated: batch,dist,insert,constraint,query,concurrent,memory,cgal\n"
        "                     or \"all\" (default; cgal needs -DQDT_BENCH_CGAL at build time)\n"
        "  --max-points N     skip every size above N\n"
        "  --csv FILE         also write all measurements to FILE\n"
        "  --help\n"
        "Example: %s --output bench.md --note \"libqdt: MSVC /O2, static library\"\n"
        "Progress goes to stderr, the report to the --output file (or to stdout).\n",
        argv0, argv0);
}

static bool valid_scenario(const std::string& s)
{
    for (size_t i = 0; i < sizeof(kScenarioNames) / sizeof(kScenarioNames[0]); ++i)
    {
        if (s == kScenarioNames[i])
        {
            return true;
        }
    }
    return false;
}

static bool parse_scenarios(const std::string& list, std::set<std::string>& out)
{
    size_t pos = 0;
    while (pos <= list.size())
    {
        size_t comma = list.find(',', pos);
        if (comma == std::string::npos)
        {
            comma = list.size();
        }
        const std::string item = list.substr(pos, comma - pos);
        pos = comma + 1;
        if (item.empty())
        {
            continue;
        }
        if (item == "all")
        {
            for (size_t i = 0; i < sizeof(kScenarioNames) / sizeof(kScenarioNames[0]); ++i)
            {
                out.insert(kScenarioNames[i]);
            }
        }
        else if (valid_scenario(item))
        {
            out.insert(item);
        }
        else
        {
            std::fprintf(stderr, "unknown scenario: %s\n", item.c_str());
            return false;
        }
    }
    return !out.empty();
}

int main(int argc, char** argv)
{
#if defined(_WIN32)
    /* Internal: measurement child of the memory scenario (see mem_child_main). */
    if (argc >= 2 && std::strcmp(argv[1], "--mem-child") == 0)
    {
        return mem_child_main(argc, argv);
    }
#endif

    Config      cfg;
    std::string scenario_arg = "all";
    bool        scenario_explicit = false;

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--quick")
        {
            cfg.quick = true;
        }
        else if (a == "--repeat" && i + 1 < argc)
        {
            cfg.reps = static_cast<size_t>(std::strtoull(argv[++i], NULL, 10));
            if (cfg.reps == 0)
            {
                std::fprintf(stderr, "--repeat needs a positive number\n");
                return 2;
            }
        }
        else if (a == "--seed" && i + 1 < argc)
        {
            cfg.seed = static_cast<uint64_t>(std::strtoull(argv[++i], NULL, 10));
        }
        else if (a == "--max-points" && i + 1 < argc)
        {
            cfg.max_points = static_cast<size_t>(std::strtoull(argv[++i], NULL, 10));
        }
        else if (a == "--scenario" && i + 1 < argc)
        {
            scenario_arg = argv[++i];
            scenario_explicit = true;
        }
        else if (a == "--csv" && i + 1 < argc)
        {
            cfg.csv_path = argv[++i];
        }
        else if (a == "--output" && i + 1 < argc)
        {
            cfg.output_path = argv[++i];
        }
        else if (a == "--note" && i + 1 < argc)
        {
            cfg.notes.push_back(argv[++i]);
        }
        else if (a == "--help" || a == "-h")
        {
            usage(argv[0]);
            return 0;
        }
        else
        {
            usage(argv[0]);
            return 2;
        }
    }
    if (!parse_scenarios(scenario_arg, cfg.scenarios))
    {
        usage(argv[0]);
        return 2;
    }

#if !defined(QDT_BENCH_CGAL)
    if (cfg.scenarios.count("cgal") != 0)
    {
        if (scenario_explicit && scenario_arg != "all")
        {
            std::fprintf(stderr, "note: built without -DQDT_BENCH_CGAL, the cgal scenario is skipped\n");
        }
        cfg.scenarios.erase("cgal");
    }
#else
    (void)scenario_explicit;
#endif
    if (cfg.scenarios.empty())
    {
        std::fprintf(stderr, "nothing to run\n");
        return 2;
    }

    build_plan(cfg);

    /* Binary mode: the report always uses LF line ends, also on Windows. */
    if (!cfg.output_path.empty())
    {
        if (!std::freopen(cfg.output_path.c_str(), "wb", stdout))
        {
            std::fprintf(stderr, "could not open %s for writing\n", cfg.output_path.c_str());
            return 2;
        }
    }

    try
    {
        print_environment(cfg);

        if (cfg.scenarios.count("batch")) bench_batch(cfg);
        if (cfg.scenarios.count("dist")) bench_dist(cfg);
        if (cfg.scenarios.count("insert")) bench_insert(cfg);
        if (cfg.scenarios.count("constraint")) bench_constraint(cfg);
        if (cfg.scenarios.count("query")) bench_query(cfg);
        if (cfg.scenarios.count("concurrent")) bench_concurrent(cfg);
        if (cfg.scenarios.count("memory")) bench_memory(cfg);
#if defined(QDT_BENCH_CGAL)
        if (cfg.scenarios.count("cgal")) bench_cgal(cfg);
#endif
    }
    catch (const std::exception& ex)
    {
        std::fprintf(stderr, "benchmark aborted: %s\n", ex.what());
        return 1;
    }

    /* A clean run adds nothing to the report; failed checks are flagged in it. */
    if (!g_warnings.empty())
    {
        std::printf("> **Warning:** some sanity checks failed. Do not publish these numbers before the cause is "
            "understood.\n>\n");
        for (size_t i = 0; i < g_warnings.size(); ++i)
        {
            std::printf("> - %s\n", g_warnings[i].c_str());
        }
        std::printf("\n");
    }
    else
    {
        std::fprintf(stderr, "  [bench] all sanity checks passed\n");
    }
    std::fflush(stdout);

    if (!cfg.csv_path.empty())
    {
        if (!write_csv(cfg.csv_path))
        {
            std::fprintf(stderr, "could not write %s\n", cfg.csv_path.c_str());
            return 1;
        }
        std::fprintf(stderr, "  [bench] wrote %zu measurements to %s\n", g_csv.size(), cfg.csv_path.c_str());
    }
    if (!cfg.output_path.empty())
    {
        std::fprintf(stderr, "  [bench] report written to %s\n", cfg.output_path.c_str());
    }
    return g_warnings.empty() ? 0 : 1;
}