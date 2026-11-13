# QDT - Quick 2D Delaunay Triangulation

QDT is a small library for fast Delaunay triangulation of 2D point sets. It is
implemented in C++14 and exposes a plain C API built around an opaque handle,
so it can be used from C, C++, and any language that can call C functions. It
supports single-threaded and parallel batch triangulation, incremental point
insertion, and constrained Delaunay triangulation (forcing given segments into
the triangulation).

## Features

- Delaunay triangulation of 2D points with `double` coordinates.
- Simple C API with an opaque handle; no dependencies beyond the C++ standard
  library.
- Robust geometric predicates: orientation and in-circle tests use a
  floating-point filter with an exact-arithmetic fallback, so rounding errors
  do not corrupt the result on degenerate or nearly degenerate input
  (collinear points, cocircular points, duplicates).
- Single-threaded and parallel batch modes.
- Incremental insertion of individual points into an existing triangulation.
- Constrained Delaunay triangulation: force segments between existing points
  to be edges, one at a time or in a batch.
- Suitable for point clouds, mesh generation, interpolation, and related
  geometry processing tasks.
- MIT licensed.

## Requirements

Using the library:

- Any C99 (or later) or C++ compiler for code that includes `qdt.h`. The header
  only depends on `<stddef.h>`.
- On Windows, see [Windows notes](#windows-notes) for the `QDT_DLL` macro.

Building the library from source:

- A C++14 (or later) compiler and a threads library (`-pthread` on Linux). The
  sources compile with GCC 13 using `-std=c++14`; C++11 is not sufficient.
- Source files: `qdt.cpp`, `GuibasStolfiDwyer_Wrapper.cpp`,
  `GuibasStolfiDwyer_Threads.cpp`, plus the headers `qdt.h` and `wrapper.hxx`.

## Building

Linux and macOS, shared library:

```sh
g++ -std=c++14 -O2 -pthread -fPIC -shared \
    qdt.cpp GuibasStolfiDwyer_Wrapper.cpp GuibasStolfiDwyer_Threads.cpp \
    -o libqdt.so
```

Static library:

```sh
g++ -std=c++14 -O2 -pthread -fPIC -c \
    qdt.cpp GuibasStolfiDwyer_Wrapper.cpp GuibasStolfiDwyer_Threads.cpp
ar rcs libqdt.a qdt.o GuibasStolfiDwyer_Wrapper.o GuibasStolfiDwyer_Threads.o
```

Linking a C program against the static library needs the C++ runtime, the math
library (the exact predicates use `fma`), and threads. Linking with `g++`
instead of `gcc` adds the first two automatically:

```sh
gcc -std=c99 example.c -I. -L. -lqdt -lstdc++ -lm -pthread -o example
```

Add the directory containing `qdt.h` to your include path in your own build.

Do not build with `-ffast-math` (GCC/Clang) or `/fp:fast` (MSVC). The exact
predicates rely on strict IEEE-754 arithmetic. On 32-bit x86, use SSE2 math
rather than x87.

### Windows notes

`qdt.h` marks every function `__declspec(dllimport)` unless `QDT_DLL` is
defined, in which case they are marked `__declspec(dllexport)`.

- When building the DLL, pass `/DQDT_DLL` to the compiler for all QDT source
  files. Defining it inside your own source file has no effect on the library.
- Code that links against the DLL includes `qdt.h` without `QDT_DLL`.
- If you build QDT as a static library or compile the sources directly into
  your application, define `QDT_DLL` for the QDT sources and for the code that
  includes `qdt.h`. Otherwise the functions are declared `dllimport` and
  linking fails.

```bat
cl /O2 /EHsc /LD /DQDT_DLL qdt.cpp GuibasStolfiDwyer_Wrapper.cpp GuibasStolfiDwyer_Threads.cpp /Fe:qdt.dll
```

## API

Include the header:

```c
#include "qdt.h"
```

### Types

```c
typedef struct qdt_vec2d
{
    double x;
    double y;
} qdt_vec2d;

typedef struct qdt_edge /* A constraint: two point indices. */
{
    size_t v0;
    size_t v1;
} qdt_edge;

typedef struct qdt_t qdt_t; /* Opaque triangulation context. */
```

### Functions

| Function | Description |
| --- | --- |
| `qdt_t* qdt_create(void)` | Creates a new, empty triangulation context. See [Error handling](#error-handling). |
| `void qdt_execute(qdt_t* ptr, const qdt_vec2d* pts, size_t n)` | Computes the Delaunay triangulation of `n` points in one batch. Replaces everything previously stored in the context, including inserted points and constraints. The points are copied; `pts` is not modified and need not stay valid after the call. |
| `void qdt_execute_parallel(qdt_t* ptr, const qdt_vec2d* pts, size_t n)` | Same as `qdt_execute`, but may use multiple threads. See [Parallel mode](#parallel-mode). |
| `size_t qdt_insert_point(qdt_t* ptr, qdt_vec2d pt)` | Appends one point to the context and updates the triangulation incrementally. Returns the total number of points in the context after the call, so the new point's index is the return value minus 1. Returns `0` if `ptr` is `NULL`. |
| `size_t qdt_insert_point_parallel(qdt_t* ptr, qdt_vec2d pt)` | Same result as `qdt_insert_point`. Only a one-time setup step uses multiple threads. See [Parallel mode](#parallel-mode). |
| `size_t qdt_insert_constraint(qdt_t* ptr, size_t v0, size_t v1)` | Forces the segment between the points with indices `v0` and `v1` to be an edge of the triangulation. Returns `1` if the edge is now a constrained edge (including when it already was an edge), `0` if the constraint was rejected or `ptr` is `NULL`. See [Constrained triangulation](#constrained-triangulation). |
| `size_t qdt_insert_constraints_batch(qdt_t* ptr, const qdt_edge* edges, size_t num_edges)` | Same as calling `qdt_insert_constraint` for each entry in order, but much faster for more than a few constraints. A rejected entry does not stop the others. Returns the number of accepted entries. |
| `size_t qdt_triangles_size(const qdt_t* ptr)` | Returns the current number of triangles. |
| `size_t qdt_get_triangle(const qdt_t* ptr, size_t T[3], size_t id)` | Writes the three point indices of triangle `id` to `T`. Returns `1` when the triangle exists, `0` otherwise (`T` is left untouched). |
| `size_t qdt_points_size(const qdt_t* ptr)` | Returns the number of points in the context: the points of the last batch call plus the points inserted since, duplicates included. Valid point indices are `0 .. qdt_points_size() - 1`. |
| `size_t qdt_get_point(const qdt_t* ptr, qdt_vec2d* pt, size_t id)` | Writes the coordinates of point `id` to `pt`. Returns `1` when the point exists, `0` otherwise (`pt` is left untouched). |
| `size_t qdt_point_is_used(const qdt_t* ptr, size_t id)` | Returns `1` if at least one triangle references point `id`. Returns `0` for a duplicate of a point with a lower index, for every point when the context has no triangles, and for an invalid `id`. |
| `void qdt_destroy(qdt_t* ptr)` | Destroys the context and releases its resources. `NULL` is allowed. |

Passing a `NULL` context to any function is safe: the functions do nothing and
return `0`.

### Results

- **Point indices.** For `qdt_execute` and `qdt_execute_parallel`, indices
  refer to the input array order. Points added with `qdt_insert_point` or
  `qdt_insert_point_parallel` are appended in insertion order. Indices stay
  valid until the next batch call.
- **Triangle ids.** Valid ids are `0 .. qdt_triangles_size() - 1`. Ids are
  **not stable**: any insertion or constraint can remove triangles and reuse
  their ids, so do not cache ids across `qdt_execute*`, `qdt_insert_point*`, or
  constraint calls.
- **Orientation.** The vertices of each triangle are listed counter-clockwise
  (positive signed area in a right-handed coordinate system). Which vertex
  comes first is unspecified, as is the order of the triangles.
- **Duplicate points.** If several points have identical coordinates, only the
  one with the lowest index appears in the triangles. The others still occupy
  an index (and count towards the value returned by `qdt_insert_point`), but
  are never referenced by a triangle.
- **Degenerate input.** With fewer than three distinct points, or with all
  points collinear, the result has zero triangles.
- **Input requirements.** Coordinates must be finite. A batch call
  (`qdt_execute*`) with a `NaN` or infinite coordinate, a `NULL` array with
  `n > 0`, or more than 100,000,000 points is rejected: the context is emptied
  (no points, no triangles) and nothing else is reported. `qdt_insert_point*`
  ignores a point with a non-finite coordinate, or one that would take the
  context beyond 100,000,000 points, and returns the unchanged number of
  points, so a return value that did not grow signals the rejection. The
  geometric tests are exact for any overall scale of the
  input and for sets that mix very different magnitudes, with one exception:
  nonzero coordinates with an absolute value below about 2.2e-308 (subnormal
  numbers) have too few significant bits, and the output may then fail to be
  Delaunay without any error being reported.

### Batch and incremental usage

QDT supports three main workflows:

1. **Batch triangulation:** call `qdt_execute` or `qdt_execute_parallel` with
   the full point array.
2. **Incremental insertion:** build the triangulation for most points with
   `qdt_execute` or `qdt_execute_parallel`, then add the remaining points one
   by one with `qdt_insert_point` or `qdt_insert_point_parallel`.
3. **Constrained triangulation:** after a batch call, force segments into the
   triangulation with `qdt_insert_constraints_batch` (or
   `qdt_insert_constraint`). See [Constrained
   triangulation](#constrained-triangulation).

Incremental insertion is meant for adding a small number of points to an
existing triangulation. Keep these costs in mind:

- The **first insertion** after a batch call builds an internal half-edge mesh.
  This takes time proportional to the size of the triangulation (about 0.7 s
  for 2.6 million points in the benchmark below, about 0.44 s with the parallel
  variant) and the mesh stays in memory until the next batch call or
  `qdt_destroy`. The first constraint call after a batch call builds the same
  mesh.
- **Every insertion** costs time proportional to the current number of points,
  not to the size of the affected neighborhood. It is not a constant-time or
  logarithmic operation. In the benchmark below it takes about 0.26 ms at
  100,000 points, 3.2 ms at 1 million points, and 8.3 ms at 2.6 million
  points.
- In the benchmark below, a full rebuild with `qdt_execute` is faster than
  inserting points one by one once more than roughly 140 points (roughly 50 to
  200 with the parallel variants, depending on the size) have to be added. If
  the point set changes by more than that, call `qdt_execute` or
  `qdt_execute_parallel` again with the complete point array.
- Inserting into a context that has no triangles yet (fewer than three
  distinct points, or all points collinear) works, but each insertion then
  triggers a full rebuild.

### Constrained triangulation

`qdt_insert_constraint(ptr, v0, v1)` forces the segment between the points
with indices `v0` and `v1` to be an edge of the triangulation. The result is
the constrained Delaunay triangulation of all points and all constraints
inserted so far: constrained edges are part of the output, no other edge
crosses them, and all other edges stay locally Delaunay. With points in
general position that triangulation is unique, so the order of the constraints
does not matter; with cocircular points, equally valid results can differ.

- A constraint is **rejected**, and the triangulation is left exactly as it
  was, when it would cross an edge that is already constrained, when the
  segment passes exactly through a third point (insert the pieces between
  consecutive points instead), when `v0` and `v1` are the same point or an
  index is not a point of the context, or when the context has no triangles.
- The order of `v0` and `v1` does not matter. The index of a duplicate point
  acts as the lowest index with the same coordinates.
- Constraints are **kept** by `qdt_insert_point*`: a point inserted exactly on
  a constrained edge splits it, and both halves stay constrained. They are
  **discarded** by the next `qdt_execute` or `qdt_execute_parallel`.
- A single `qdt_insert_constraint` call finds its two end points in a table
  that is built together with the internal mesh and kept up to date by every
  insertion, so the lookup takes constant time and does not grow with the
  number of points. `qdt_insert_constraints_batch` gives the same result and is
  the convenient way to insert many constraints in one call. The cost of a
  constraint grows with the number of triangles its segment crosses. (The
  benchmark below was measured before the lookup table existed, when a single
  call searched linearly; its single-call timings no longer apply and have not
  been re-measured.)

```c
/* A flat quadrilateral: the Delaunay triangulation uses the diagonal 1-3.
   The constraint 0-2 forces the other diagonal. */
qdt_vec2d pts[] = { {0.0, 0.0}, {4.0, 0.0}, {3.0, 1.0}, {1.0, 1.0} };

qdt_t* q = qdt_create();
qdt_execute(q, pts, 4);          /* triangles: (3, 0, 1) and (1, 2, 3) */
qdt_insert_constraint(q, 0, 2);  /* returns 1; triangles: (2, 0, 1) and (0, 2, 3) */
qdt_destroy(q);
```

### Parallel mode

`qdt_execute_parallel` produces the same set of triangles as `qdt_execute`,
but not necessarily in the same order.

- The thread count is `std::thread::hardware_concurrency()`, capped at 64, and
  is not configurable. Threads exist only for the duration of the call.
- With fewer than 200,000 distinct points (or a single hardware thread), the
  call runs on the calling thread.
- In the current implementation, a subproblem is split between threads only if
  it has at least 524,288 points, and at most eight subproblems are built
  concurrently. One million points are therefore processed as two parallel
  subproblems, and about 2.1 million points or more as eight, no matter how
  many hardware threads exist. Inputs with fewer than 524,288 distinct points
  get no parallel speedup at all (this includes 1,000,000 points of which half
  are duplicates).
- For inputs of up to a few hundred thousand points, expect no benefit from
  the parallel version: it uses a separate implementation that pre-allocates
  working storage for at least 65,536 edge records (about 1.5 MiB) even for
  tiny inputs. In the benchmark below it takes about 9% longer than
  `qdt_execute` at 100,000 points, and much longer in relative terms for tiny
  inputs. Use `qdt_execute` for such inputs.
- If the parallel implementation cannot run (its pre-allocated edge storage is
  exhausted, or the worker threads cannot be started), the call falls back to
  the serial algorithm and still returns the complete result.

`qdt_insert_point_parallel` gives the same result as `qdt_insert_point`. The
insertion itself is single-threaded; only the one-time mesh build performed by
the first insertion after a batch call is parallel (and only when the
triangulation has at least 50,000 triangles). This matches the benchmark
below: `insert(0)` is faster in the parallel run, while later insertions take
about the same time in both.

### Error handling

QDT has no error codes and does not catch C++ exceptions at the API boundary.
Allocation failures surface as a `std::bad_alloc` exception thrown from
`qdt_create`, `qdt_execute*`, or `qdt_insert_point*`. `qdt_create` does not
return `NULL`. Letting an exception propagate into C code is undefined
behavior, so treat out-of-memory as fatal when calling from C.

### Thread safety

A single context must not be used from several threads at the same time.
Different contexts share no state and can be used concurrently.

### Resource use

Measured with the benchmark below (Windows, committed private bytes per input
point, 1 to 2.6 million points):

- A context holds about 72 bytes per input point after a batch call: 16 for
  the point copy, 8 for the index map, and 48 for the triangles (two per
  point). Small inputs show somewhat more, for example 96 bytes per point at
  100,000 points.
- At its peak, a batch call uses about 150 bytes per point above the level
  before the call (about 200 with `qdt_execute_parallel`); about 72 of these
  remain afterwards.
- The mesh built by the first insertion or constraint call stays in memory and
  adds about 400 bytes per point (committed memory; the point-to-vertex lookup
  table added later adds 4 more bytes per point, not included in the figures
  here). While it is being built,
  the usage peaks at about 515 bytes per point above the level after the batch
  call. The committed figure counts storage that is reserved but not fully
  used. On Linux, where only touched pages count, the resident set grew by
  about 180 bytes per point for the mesh in a measurement with 100,000 points.
  Call `qdt_execute*` again or `qdt_destroy` to release the mesh.
- Internal indices are 32-bit, so inputs of more than 100,000,000 points are
  rejected (see [Results](#results)).

## Example

The complete, compilable program is [`example/example.c`](example/example.c).
These are the calls it makes, in order (`pts` is an array of `n` points):

```c
qdt_t* q = qdt_create();

/* Batch-triangulate most of the points. */
qdt_execute(q, pts, n - 2);

/* Insert the remaining points incrementally. */
qdt_insert_point(q, pts[n - 2]);
qdt_insert_point(q, pts[n - 1]);

/* Use qdt_execute_parallel(q, pts, n) for a full parallel batch. */
/* Use qdt_insert_point_parallel(q, pt) for parallel insertion. */

for (size_t i = 0; i < qdt_triangles_size(q); ++i)
{
    size_t tri[3];
    if (qdt_get_triangle(q, tri, i))
    {
        /* tri[0], tri[1], tri[2] index into pts, counter-clockwise. */
    }
}

qdt_destroy(q);
```

The program triangulates six points. Points 0-3 form a square (all four are
cocircular), point 4 is its center, and point 5 lies exactly on the edge between
points 4 and 2.

Example output (six points, four of them on the convex hull, so
`2 * 6 - 2 - 4 = 6` triangles):

```text
Triangles: 6
Triangle 0: 1, 3, 4
Triangle 1: 0, 1, 4
Triangle 2: 2, 0, 5
Triangle 3: 0, 4, 5
Triangle 4: 4, 3, 5
Triangle 5: 3, 2, 5
```

The order of the triangles, and which vertex of each triangle comes first, may
differ between versions and platforms.

## Performance

The figures below were produced by `test/qdt_bench.cpp`. Timings are wall-clock milliseconds around the library calls only (point generation and the creation and destruction of contexts are not included); tables show medians. Absolute numbers depend on the CPU, compiler, compiler flags and thread count, so the relative comparisons are more useful than the absolute values.

### Test environment

- Date: 2026-11-10 01:02 UTC
- CPU: Intel(R) Core(TM) i7-8700 CPU @ 3.20GHz
- Cores / hardware threads: 6 / 12 (`qdt_execute_parallel` uses min(threads, 64))
- Memory: 31.9 GiB
- OS: Windows 11 25H2 (build 26200), x64
- Compiler (benchmark program): MSVC 19.51.36260, 64-bit
- Parameters: repeat = 5, seed = 20260607

### Batch triangulation: scaling with the number of points

Uniform random points; medians of 5 runs (3 for the largest sizes). `qdt_execute_parallel` uses several threads only from 200,000 distinct points on; below that it runs on the calling thread. Inputs below 50,000 points are timed as the average of several consecutive calls on one context.

| Points    | Triangles | Serial (ms) | Parallel (ms) | Speedup | Serial (Mpts/s) | Parallel (Mpts/s) |
| --------- | --------: | ----------: | ------------: | ------: | --------------: | ----------------: |
| 100       |       185 |      0.0284 |        0.4826 |   0.06x |            3.52 |              0.21 |
| 1,000     |     1,976 |      0.5470 |         1.030 |   0.53x |            1.83 |              0.97 |
| 10,000    |    19,968 |       5.914 |         6.868 |   0.86x |            1.69 |              1.46 |
| 100,000   |   199,958 |       63.38 |         69.14 |   0.92x |            1.58 |              1.45 |
| 1,000,000 | 1,999,963 |       702.5 |         522.9 |   1.34x |            1.42 |              1.91 |
| 2,624,409 | 5,248,777 |      1951.9 |         913.3 |   2.14x |            1.34 |              2.87 |

### Batch triangulation: point distributions

1,000,000 points per case; medians of 5 runs. Triangle counts differ between cases because the number of convex hull points and of duplicate points differs.

| Distribution             | Triangles | Serial (ms) | Parallel (ms) | Speedup | Serial (Mpts/s) |
| ------------------------ | --------: | ----------: | ------------: | ------: | --------------: |
| uniform random           | 1,999,953 |       708.1 |         502.6 |   1.41x |            1.41 |
| uniform, presorted by x  | 1,999,959 |       623.9 |         422.1 |   1.48x |            1.60 |
| gaussian                 | 1,999,980 |       726.3 |         518.2 |   1.40x |            1.38 |
| clustered (20 gaussians) | 1,999,977 |       719.8 |         507.6 |   1.42x |            1.39 |
| regular grid (shuffled)  | 1,996,002 |      1131.4 |         721.5 |   1.57x |            0.88 |
| jittered grid (shuffled) | 1,999,961 |       690.3 |         493.1 |   1.40x |            1.45 |
| all points on a circle   | 1,000,017 |       556.5 |         404.9 |   1.37x |            1.80 |
| flat band (aspect 1e6:1) | 1,999,962 |      1614.3 |        1062.1 |   1.52x |            0.62 |
| 50% duplicate points     |   999,964 |       419.0 |         448.3 |   0.93x |            2.39 |

### Incremental insertion

A base triangulation of N uniform random points, followed by 50 single point insertions. The first insertion builds the internal mesh and is much slower than the later ones. Break-even is the estimated number of inserted points at which a full rebuild becomes cheaper: 1 + (build - first insert) / mean(later inserts). Serial = `qdt_execute` + `qdt_insert_point`; parallel = `qdt_execute_parallel` + `qdt_insert_point_parallel`. Medians of up to 3 runs.

| Base points | Mode     | Build (ms) | 1st insert (ms) | Later median (ms) | Later mean (ms) | Later p95 (ms) | Break-even (points) |
| ----------- | -------- | ---------: | --------------: | ----------------: | --------------: | -------------: | ------------------: |
| 10,000      | serial   |      5.882 |           1.753 |            0.0283 |          0.0291 |         0.0334 |                ~143 |
| 10,000      | parallel |      6.769 |           1.715 |            0.0281 |          0.0289 |         0.0327 |                ~176 |
| 100,000     | serial   |      63.76 |           20.48 |            0.2631 |          0.2950 |         0.3746 |                ~148 |
| 100,000     | parallel |      69.82 |           14.73 |            0.2499 |          0.2825 |         0.3316 |                ~196 |
| 1,000,000   | serial   |      704.6 |           230.4 |             3.187 |           3.428 |          3.502 |                ~139 |
| 1,000,000   | parallel |      507.5 |           158.4 |             3.184 |           3.404 |          3.401 |                ~104 |
| 2,624,409   | serial   |     1886.6 |           672.3 |             8.309 |           8.936 |          8.964 |                ~137 |
| 2,624,409   | parallel |      905.6 |           443.0 |             8.235 |           8.806 |          8.631 |                 ~54 |

#### Per-insert detail, non-parallel, base 2,624,409 points

| Operation               | Triangles | Time (ms) |
| ----------------------- | --------: | --------: |
| Full (2,624,419 points) | 5,248,798 |    1889.4 |
| Most, 10 points left    | 5,248,778 |    1887.5 |
| insert(0)               | 5,248,780 |     696.1 |
| insert(1)               | 5,248,782 |     8.838 |
| insert(2)               | 5,248,784 |     8.753 |
| insert(3)               | 5,248,786 |     8.938 |
| insert(4)               | 5,248,788 |     8.665 |
| insert(5)               | 5,248,790 |     8.861 |
| insert(6)               | 5,248,792 |     8.706 |
| insert(7)               | 5,248,794 |     8.554 |
| insert(8)               | 5,248,796 |     8.700 |
| insert(9)               | 5,248,798 |     8.836 |

#### Per-insert detail, parallel, base 2,624,409 points

| Operation               | Triangles | Time (ms) |
| ----------------------- | --------: | --------: |
| Full (2,624,419 points) | 5,248,798 |     907.1 |
| Most, 10 points left    | 5,248,778 |     908.2 |
| insert(0)               | 5,248,780 |     442.2 |
| insert(1)               | 5,248,782 |     8.286 |
| insert(2)               | 5,248,784 |     8.256 |
| insert(3)               | 5,248,786 |     8.294 |
| insert(4)               | 5,248,788 |     8.220 |
| insert(5)               | 5,248,790 |     8.223 |
| insert(6)               | 5,248,792 |     8.294 |
| insert(7)               | 5,248,794 |     8.413 |
| insert(8)               | 5,248,796 |     8.234 |
| insert(9)               | 5,248,798 |     8.457 |

### Constrained Delaunay triangulation

A base triangulation of N uniform random points, then constraints from three sets: "short" (segments between nearby points, none of them crossing; some are already Delaunay edges), "long chords" (segments spanning the whole domain, each crossing about sqrt(N) triangles) and "ring" (a closed polygon through the points of a thin annulus). The first constraint is inserted with `qdt_insert_constraint` and includes the one-time internal mesh build; the others are inserted either with one `qdt_insert_constraints_batch` call or with a loop of `qdt_insert_constraint` calls. Every constraint of these sets is valid and must be accepted. Medians of up to 3 runs.

| Points    | Constraint set | Constraints | 1st call + mesh build (ms) | Batch (ms) | Batch (us/constraint) | Single calls (ms) | Single (us/constraint) |        Accepted |
| --------- | -------------- | ----------: | -------------------------: | ---------: | --------------------: | ----------------: | ---------------------: | --------------: |
| 10,000    | short          |         492 |                      1.475 |     0.6266 |                  1.28 |             3.039 |                   6.19 |       492 / 492 |
| 10,000    | long chords    |          35 |                      1.632 |      2.256 |                  66.3 |             2.458 |                   72.3 |         35 / 35 |
| 10,000    | ring           |       1,842 |                      1.494 |      2.054 |                  1.12 |             10.99 |                   5.97 |   1,842 / 1,842 |
| 100,000   | short          |       4,839 |                      16.41 |      7.391 |                  1.53 |             293.5 |                   60.7 |   4,839 / 4,839 |
| 100,000   | long chords    |         100 |                      17.34 |      23.92 |                   242 |             29.71 |                    300 |       100 / 100 |
| 100,000   | ring           |       1,942 |                      16.17 |      1.259 |                  0.65 |             114.8 |                   59.1 |   1,942 / 1,942 |
| 1,000,000 | short          |      20,000 |                      195.3 |      40.17 |                  2.01 |           27284.0 |                   1364 | 20,000 / 20,000 |
| 1,000,000 | long chords    |         100 |                      196.1 |      105.0 |                  1060 |             240.3 |                   2427 |       100 / 100 |
| 1,000,000 | ring           |       2,006 |                      194.9 |      5.215 |                  2.60 |            2642.5 |                   1318 |   2,006 / 2,006 |

Accepted is shown as batch / single-call counts.

### Reading the result

`qdt_get_triangle` over all triangles of a uniform random point set.

| Points    | Triangles | Time (ms) | ns per call | Mcalls/s |
| --------- | --------: | --------: | ----------: | -------: |
| 1,000,000 | 1,999,957 |     5.244 |         2.6 |    381.4 |

### Independent contexts on several threads

Every thread triangulates the same 100,000 uniform random points with `qdt_execute` on its own context (the contexts share no state). Scaling is the aggregate throughput divided by the throughput of one thread; efficiency is scaling divided by the thread count. Hardware threads: 12. Medians of 5 runs.

| Threads | Wall time (ms) | Aggregate (Mpts/s) | Scaling | Efficiency |
| ------- | -------------: | -----------------: | ------: | ---------: |
| 1       |          64.60 |               1.55 |   1.00x |       100% |
| 2       |          66.21 |               3.02 |   1.95x |        98% |
| 4       |          75.55 |               5.29 |   3.42x |        86% |
| 8       |          95.29 |               8.40 |   5.42x |        68% |
| 12      |          112.5 |              10.67 |   6.89x |        57% |

### Memory use

Committed private bytes (PrivateUsage) measured around the library calls, in a freshly started process per measurement; the input point array is part of the baseline. "Retained" is the growth after the call. "Peak" is the highest level during the call minus the level before it. Mesh = the internal half-edge mesh built by the first insertion after a batch call. Bytes are per input point.

| Points    | Mode     | Retained after batch (B/pt) | Peak during batch (B/pt) | Mesh retained (B/pt) | Peak during mesh build (B/pt) | Retained, batch + mesh (MiB) |
| --------- | -------- | --------------------------: | -----------------------: | -------------------: | ----------------------------: | ---------------------------: |
| 100,000   | serial   |                        96.4 |                    152.5 |                401.4 |                         509.6 |                         47.5 |
| 100,000   | parallel |                        96.3 |                    216.2 |                409.5 |                         517.8 |                         48.2 |
| 1,000,000 | serial   |                        72.2 |                    152.3 |                396.8 |                         517.3 |                        447.2 |
| 1,000,000 | parallel |                        72.2 |                    202.5 |                396.8 |                         521.3 |                        447.3 |
| 2,624,409 | serial   |                        72.1 |                    152.3 |                396.8 |                         513.3 |                       1173.7 |
| 2,624,409 | parallel |                        72.2 |                    201.2 |                396.8 |                         513.3 |                       1173.7 |

### Notes on the results

- **Batch.** Serial throughput falls from about 3.5 million points per second
  for tiny inputs to about 1.3 to 1.4 million for millions of points.
  2,624,409 points (5,248,777 triangles) take 1.95 s with `qdt_execute` and
  0.91 s with `qdt_execute_parallel` (2.14x). The parallel call only pays off
  for large inputs; see [Parallel mode](#parallel-mode).
- **Distributions.** Compared with uniform random points (708 ms for 1,000,000
  points), presorted input is 12% faster and points on a circle are 21%
  faster; gaussian, clustered, and jittered-grid inputs are within 3%; a
  regular grid with many cocircular points is 1.6x slower, and the very thin
  band is 2.3x slower. The input with 50% duplicates is the fastest in absolute
  terms because only 500,000 distinct points are triangulated. That is below
  the 524,288-point threshold for splitting work between threads, so it gets
  no parallel speedup (0.93x).
- **Incremental insertion.** The first insertion costs about a third of a
  serial build (672 ms against 1887 ms at 2,624,409 points) and about half of a
  parallel build. Later insertions grow roughly linearly with the size: 0.028
  ms at 10,000 points, 0.26 ms at 100,000, 3.2 ms at 1,000,000, and 8.3 ms at
  2,624,409. Building most of the points and then inserting the last 10 one by
  one takes 1.4x (serial) or 1.6x (parallel) as long as triangulating all of
  them at once, so incremental insertion only pays off for a small number of
  points.
- **Constraints.** (Measured before the end-point lookup table was added: the
  comparison of single calls with the batch call below no longer applies.) At
  1,000,000 points the first constraint call, which
  includes building the internal mesh, takes about 195 ms (28% of the 702 ms
  batch triangulation). After that, `qdt_insert_constraints_batch` inserts
  short constraints at about 2 µs each, while the same 20,000 short constraints
  as single `qdt_insert_constraint` calls take 27 s (680 times longer). The
  factor is about 40 at 100,000 points and about 5 at 10,000 points. Long
  chords, each crossing about sqrt(N) triangles, cost about 1 ms each at
  1,000,000 points; the crossing dominates there, so the batch call is only
  2.3 times faster.
- **Several contexts.** Independent contexts scale well: 3.4x on 4 threads and
  6.9x on all 12 hardware threads of this 6-core CPU (the second thread of a
  core adds less than a full core).
- **Reading the result.** `qdt_get_triangle` takes about 2.6 ns per call.

### Reproducing the benchmark

`test/qdt_bench.cpp` uses only the public API and builds on Linux and Windows.
It writes the report above as Markdown. Progress goes to stderr, so use
`--output` rather than shell redirection (Windows PowerShell 5.1 writes
redirected output as UTF-16).

```sh
# Linux, with libqdt.a built as described in "Building"
g++ -std=c++14 -O2 -DNDEBUG -pthread -Ilibqdt test/qdt_bench.cpp libqdt.a -lm -o qdt_bench
./qdt_bench --output bench.md --csv bench.csv --note "libqdt <commit>, g++ -O2, static library"
```

```bat
:: Windows, MSVC; qdt.lib built as a static library with /DQDT_DLL
cl /std:c++14 /O2 /EHsc /DNDEBUG /DQDT_DLL /Ilibqdt test\qdt_bench.cpp qdt.lib
qdt_bench --output bench.md --csv bench.csv --note "libqdt <commit>, MSVC /O2, static library"
```

Use an optimized build of the library, close other programs, and use a
high-performance power plan. A full run takes several minutes; `--quick` runs
small sizes only, and `--scenario batch,insert,...` selects scenarios. The
first-insertion time varies the most between runs. The benchmark also has an
optional comparison with CGAL (build with `-DQDT_BENCH_CGAL`).

<!-- QDT-VERIFY-BEGIN -->
## Correctness tests

The results below were produced by `test/qdt_verify.cpp`, which uses only the public C API. It runs randomized and structured test cases and checks every triangulation with exact integer arithmetic, independent of the predicates inside the library: the inputs have integer coordinates with an absolute value of at most 2^28, which `double` represents exactly. A case is one input together with the operations applied to it (for example a point set, then insertions and constraints). The context is verified at the end of the case, and in the incremental scenarios also in between.

### Test environment

- Date: 2026-11-13 02:10 UTC
- CPU: Intel(R) Core(TM) i7-8700 CPU @ 3.20GHz
- Hardware threads: 12
- Memory: 31.9 GiB
- OS: Windows 11 25H2 (build 26200), x64
- Compiler (test program): MSVC 19.51.36260, 64-bit; exact arithmetic: built-in 160-bit integers
- Parameters: iterations = 3000, seed = 0

### Results

**All 16 scenarios passed:** 32,297 cases, 0 failures, 71.85 s.

| Scenario | What is checked | README | Cases | Failures | Time (s) |
| --- | --- | --- | ---: | ---: | ---: |
| `api` | API contract: NULL context, out-of-range ids, query functions, duplicates, degenerate contexts, constraint end points, replacement by the next batch call | [Functions](#functions), [Results](#results) | 46 | 0 | 0.00 |
| `examples` | The README examples (quadrilateral with a constraint, six-point program); the input array is copied and not modified | [Example](#example), [Constrained triangulation](#constrained-triangulation) | 12 | 0 | 0.00 |
| `guards` | Invalid input: NaN, infinity, NULL array, more than 100,000,000 points | [Results](#results) | 22 | 0 | 0.00 |
| `structured` | Degenerate and structured inputs: grids, collinear points, duplicates, fewer than three points, coordinates up to 2^28 | [Results](#results) | 160 | 0 | 0.38 |
| `batch` | Random small point sets with many duplicates and cocircular points, serial and parallel batch call | [Features](#features), [Results](#results) | 6,000 | 0 | 4.37 |
| `magnitudes` | One point set that mixes magnitudes from 1 to 2^28, including clusters far from the origin | [Results](#results) | 3,000 | 0 | 2.75 |
| `scale` | Scaling the whole input by 2^-1000 ... 2^900 must not change the result | [Results](#results) | 1,200 | 0 | 0.59 |
| `circle` | Nearly cocircular points, with insertions on top | [Features](#features) | 400 | 0 | 2.78 |
| `insert` | Insertion into small triangulations, also into a context that was never filled: duplicates, points on edges, serial and parallel | [Batch and incremental usage](#batch-and-incremental-usage) | 6,000 | 0 | 15.04 |
| `insertlarge` | Insertion into triangulations of a few thousand points: duplicates, points on edges, points outside the hull | [Batch and incremental usage](#batch-and-incremental-usage) | 31 | 0 | 1.66 |
| `parmesh` | First insertion builds the internal mesh in parallel (70,000 points, at least 50,000 triangles), then insertions and constraints | [Parallel mode](#parallel-mode) | 12 | 0 | 14.75 |
| `cdt` | Constraints against an independent oracle for the rules: accepted, rejected, unchanged after a rejection; single calls and batch | [Constrained triangulation](#constrained-triangulation) | 6,000 | 0 | 0.98 |
| `cdtins` | The same with insertions in between: constraints are kept, and split by points inserted on them | [Constrained triangulation](#constrained-triangulation) | 6,000 | 0 | 3.02 |
| `cdtequiv` | A batch call equals the single calls in the same order; in general position the order of the constraints does not matter | [Constrained triangulation](#constrained-triangulation) | 1,001 | 0 | 0.26 |
| `big` | 600,000 to 700,000 points: serial and parallel batch call give valid results and the same triangles | [Parallel mode](#parallel-mode) | 3 | 0 | 22.90 |
| `contexts` | Independent contexts used from several threads at the same time | [Thread safety](#thread-safety) | 2,410 | 0 | 2.36 |

Notes:

- `examples`: The listing of the six-point example is reproduced exactly, in the same order.
- `big`: Distributions: uniform, a 1500 x 1500 grid with many duplicates (counts compared, not sets: cocircular ties), a flat band.
- `contexts`: 8 threads, each with its own contexts (batch, parallel batch and insertions). Plus two parallel batch calls at the same time on 600,000 points.

### What is checked in every case

- Every triangle has three distinct vertices in strictly counter-clockwise order and refers only to the lowest index of a group of duplicate points.
- Every directed edge occurs once, so triangles do not overlap; the areas of the triangles add up to the area of the convex hull of the points; no hull vertex lies to the right of a boundary edge.
- Every interior edge is locally Delaunay, except for constrained edges. Without constraints, small inputs are also checked against the empty-circle property for every point.
- Every accepted constraint is an edge of the triangulation, and a rejected constraint leaves the triangulation unchanged.
- `qdt_points_size`, `qdt_get_point` and `qdt_point_is_used` agree with the points that were given and with the triangles.

### Not covered

- The fallback of `qdt_execute_parallel` to the serial algorithm (exhausted edge storage, threads that cannot be started) and out-of-memory behavior: these paths cannot be triggered through the public API.
- Inputs of more than 700,000 points, and the speed and memory figures (see [Performance](#performance)).
- Coordinates outside the range of the exact check: integers up to 2^28, scaled as a whole by powers of two. Subnormal coordinates are documented as unsupported.
- Compiling `qdt.h` as C, and every compiler and platform other than the one listed above.
- Using one context from several threads, which is not supported.

### Reproducing the tests

`test/qdt_verify.cpp` uses only the public API. `--readme README.md` replaces this section (the text between the two `QDT-VERIFY` marker comments, inserted before "Current limitations" the first time) and is refused when a check fails. The program exits with status 1 if a check fails. `--quick` runs a short smoke test and `--scenario LIST` selects scenarios; neither updates the README unless `--allow-partial` is given.

```sh
# Linux, with libqdt.a built as described in "Building"
g++ -std=c++14 -O2 -DNDEBUG -pthread -Ilibqdt test/qdt_verify.cpp libqdt.a -lm -o qdt_verify
./qdt_verify --readme README.md --note "libqdt <commit>, g++ -O2, static library"
```

```bat
:: Windows, MSVC; qdt.lib built as a static library with /DQDT_DLL
cl /std:c++14 /O2 /EHsc /DNDEBUG /DQDT_DLL /Ilibqdt test\qdt_verify.cpp qdt.lib
qdt_verify --readme README.md --note "libqdt <commit>, MSVC /O2, static library"
```
<!-- QDT-VERIFY-END -->

## Current limitations

- Incremental insertion is intended for adding a small number of points to an
  existing triangulation. Its per-point cost is proportional to the current
  number of points. If the point set changes substantially, call
  `qdt_execute` or `qdt_execute_parallel` again with the complete updated
  point array.
- The current API does not provide deletion of points or direct editing of
  existing vertices, and constraints cannot be removed individually. They are
  discarded only by the next `qdt_execute` or `qdt_execute_parallel`.
- A constraint that crosses an already constrained edge, or passes exactly
  through a third point, is rejected (see [Constrained
  triangulation](#constrained-triangulation)). Constraints join existing
  points; they never add points. The cost of one `qdt_insert_constraint` call
  depends on the number of triangles its segment crosses, not on the number of
  points.
- Triangle ids are not stable across insertions or constraints, and the order
  of triangles and of the vertices within a triangle is unspecified. Point
  indices are stable; keep the input point order and record insertion order if
  you need to map results back to your data.
- Parallel functions may use multiple threads, but actual speedup depends on
  hardware, compiler, thread count, and workload size. In the current
  implementation, only inputs of at least 524,288 distinct points get a speedup
  and at most eight subproblems run concurrently (see [Parallel
  mode](#parallel-mode)). The thread count cannot be configured.
- Errors are reported as C++ exceptions, not return codes (see
  [Error handling](#error-handling)).
- Inputs must contain only finite coordinates, and a context must not be shared
  between threads without external synchronization.

## License

QDT is released under the MIT License. See [LICENSE](LICENSE).