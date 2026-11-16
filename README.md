# QDT - Quick 2D Delaunay Triangulation

QDT is a small library for fast Delaunay triangulation of 2D point sets. It is
implemented in C++14 and exposes a plain C API built around an opaque handle,
so it can be used from C, C++, and any language that can call C functions. It
supports single-threaded and parallel batch triangulation, incremental point
insertion, and constrained Delaunay triangulation (forcing given segments into
the triangulation). Every function that can fail reports the reason as a status
code.

## Features

- Delaunay triangulation of 2D points with `double` coordinates.
- Simple C API with an opaque handle; no dependencies beyond the C++ standard
  library.
- Errors are reported as `qdt_status` return codes. No C++ exception ever
  crosses the C API: out-of-memory and other internal errors are handled inside
  the library (see [Error handling](#error-handling)), so every function is
  safe to call from C.
- Robust geometric predicates: orientation and in-circle tests use a
  floating-point filter with an exact-arithmetic fallback, so rounding errors
  do not corrupt the result on degenerate or nearly degenerate input
  (collinear points, cocircular points, duplicates).
- Single-threaded and parallel batch modes.
- Speed on par with CGAL: in the benchmark below, serial batch triangulation
  takes about as long as CGAL's `Delaunay_triangulation_2` (QDT is within 6% of
  CGAL from 10,000 to 2.6 million points), and the parallel variant is 1.35x
  faster than CGAL at 1 million points and 2.0x faster at 2.6 million points
  (see [Comparison with CGAL](#comparison-with-cgal)).
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
  C++ exceptions must be enabled (the default): the library throws and catches
  them internally, so do not build it with `-fno-exceptions` or without `/EHsc`.
- Source files: `qdt.cpp`, `GuibasStolfiDwyer_Wrapper.cpp`,
  `GuibasStolfiDwyer_Threads.cpp`, plus the headers `qdt.h`, `wrapper.hxx`,
  `builder.hxx` and `robust.hxx`. Only `qdt.h` is needed by code that uses the
  library.

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

### CMake

The repository has a `CMakeLists.txt` that builds the library from the sources
in `libqdt/`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix /path/to/prefix
```

- `BUILD_SHARED_LIBS` selects a shared library; the default is a static one.
- `QDT_BUILD_BENCHMARK=ON` also builds the `qdt_bench` executable. This target
  does not include the optional CGAL comparison; build `qdt_bench` by hand for
  that (see [Reproducing the benchmark](#reproducing-the-benchmark)).
- The installed package exports the target `QDT::qdt`:
  `find_package(qdt CONFIG REQUIRED)` and
  `target_link_libraries(app PRIVATE QDT::qdt)`. A project that only uses the
  C language must also enable C++ (`project(app C CXX)`) so that the static
  library is linked with the C++ runtime.
- On Windows the build defines `QDT_DLL` as described below.

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

typedef enum qdt_status /* Result of a function that can fail; see below. */
{
    QDT_OK = 0,
    QDT_ERR_NULL_CONTEXT = -1,
    /* ... the other codes are negative, see "Status codes" */
} qdt_status;

typedef struct qdt_t qdt_t; /* Opaque triangulation context. */
```

### Functions

| Function | Description |
| --- | --- |
| `qdt_t* qdt_create(void)` | Creates a new, empty triangulation context. Returns `NULL` if the memory for it cannot be allocated. |
| `qdt_status qdt_execute(qdt_t* ptr, const qdt_vec2d* pts, size_t n)` | Computes the Delaunay triangulation of `n` points in one batch. Replaces everything previously stored in the context, including inserted points and constraints. The points are copied; `pts` is not modified and need not stay valid after the call. Invalid input and running out of memory return an error code and leave the context empty (see [Error handling](#error-handling)). |
| `qdt_status qdt_execute_parallel(qdt_t* ptr, const qdt_vec2d* pts, size_t n)` | Same as `qdt_execute`, but may use multiple threads. See [Parallel mode](#parallel-mode). |
| `qdt_status qdt_insert_point(qdt_t* ptr, qdt_vec2d pt, size_t* index)` | Appends one point to the context and updates the triangulation incrementally. `index` is optional (may be `NULL`); on success it receives the index of the new point, duplicates included. It is not written when the call fails. A point with a non-finite coordinate, or one that would take the context beyond 100,000,000 points, is rejected and nothing changes. |
| `qdt_status qdt_insert_point_parallel(qdt_t* ptr, qdt_vec2d pt, size_t* index)` | Same result as `qdt_insert_point`. Only a one-time setup step uses multiple threads. See [Parallel mode](#parallel-mode). |
| `qdt_status qdt_insert_constraint(qdt_t* ptr, size_t v0, size_t v1)` | Forces the segment between the points with indices `v0` and `v1` to be an edge of the triangulation. Returns `QDT_OK` if the edge is now a constrained edge (including when it already was an edge), or the reason why the constraint was rejected. See [Constrained triangulation](#constrained-triangulation). |
| `qdt_status qdt_insert_constraints_batch(qdt_t* ptr, const qdt_edge* edges, size_t num_edges, qdt_status* results)` | Same as calling `qdt_insert_constraint` for each entry in order. A rejected entry does not stop the others. `results` is optional (may be `NULL`); if given, it must have `num_edges` entries and entry `i` receives the status of `edges[i]`. Returns `QDT_OK` if every entry was accepted (also for `num_edges == 0`), otherwise the status of the first rejected entry. |
| `size_t qdt_triangles_size(const qdt_t* ptr)` | Returns the current number of triangles. |
| `size_t qdt_get_triangle(const qdt_t* ptr, size_t T[3], size_t id)` | Writes the three point indices of triangle `id` to `T`. Returns `1` when the triangle exists, `0` otherwise (`T` is left untouched). |
| `size_t qdt_points_size(const qdt_t* ptr)` | Returns the number of points in the context: the points of the last batch call plus the points inserted since, duplicates included. Valid point indices are `0 .. qdt_points_size() - 1`. |
| `size_t qdt_get_point(const qdt_t* ptr, qdt_vec2d* pt, size_t id)` | Writes the coordinates of point `id` to `pt`. Returns `1` when the point exists, `0` otherwise (`pt` is left untouched). |
| `size_t qdt_point_is_used(const qdt_t* ptr, size_t id)` | Returns `1` if at least one triangle references point `id`. Returns `0` for a duplicate of a point with a lower index, for every point when the context has no triangles, and for an invalid `id`. |
| `const char* qdt_status_string(qdt_status status)` | Returns a short English description of a status code. The string has static storage duration and is never `NULL`; an unknown value gives `"unknown status"`. |
| `void qdt_destroy(qdt_t* ptr)` | Destroys the context and releases its resources. `NULL` is allowed. |

Passing a `NULL` context to any function is safe. The functions that return a
`qdt_status` return `QDT_ERR_NULL_CONTEXT` and change nothing, the query
functions return `0`, and `qdt_destroy` does nothing. Unlike the functions that
return a `qdt_status`, the query functions return `1` for success.

### Status codes

`QDT_OK` is `0` and every error is negative, so `status != QDT_OK` (or
`status < 0`) tests for failure.

| Code | Value | Meaning |
| --- | ---: | --- |
| `QDT_OK` | 0 | Success. |
| `QDT_ERR_NULL_CONTEXT` | -1 | The context pointer is `NULL`. |
| `QDT_ERR_NULL_ARGUMENT` | -2 | A required pointer argument is `NULL`: `pts` with `n > 0`, or `edges` with `num_edges > 0`. |
| `QDT_ERR_TOO_MANY_POINTS` | -3 | More than 100,000,000 points. |
| `QDT_ERR_NOT_FINITE` | -4 | A coordinate is NaN or infinite. |
| `QDT_ERR_INVALID_INDEX` | -5 | A constraint end point is not a point of the context. |
| `QDT_ERR_SAME_POINT` | -6 | The end points of a constraint are the same point (the same index, or points with identical coordinates). |
| `QDT_ERR_NO_TRIANGLES` | -7 | A constraint was requested, but the context has no triangles (fewer than three distinct points, or all points collinear). |
| `QDT_ERR_CROSSES_CONSTRAINT` | -8 | The constraint would cross an edge that is already constrained. |
| `QDT_ERR_THROUGH_POINT` | -9 | The constraint segment passes exactly through a third point. |
| `QDT_ERR_OUT_OF_MEMORY` | -10 | Memory ran out. The context is emptied. |
| `QDT_ERR_INTERNAL` | -11 | Unexpected internal failure. The context is emptied. |

What an error does to the context is described under
[Error handling](#error-handling).

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
  an index (and count towards `qdt_points_size`), but are never referenced by a
  triangle.
- **Degenerate input.** With fewer than three distinct points, or with all
  points collinear, the result has zero triangles.
- **Input requirements.** Coordinates must be finite. A batch call
  (`qdt_execute*`) is rejected, and the context emptied (no points, no
  triangles), when the input has more than 100,000,000 points
  (`QDT_ERR_TOO_MANY_POINTS`), a `NULL` array with `n > 0`
  (`QDT_ERR_NULL_ARGUMENT`), or a `NaN` or infinite coordinate
  (`QDT_ERR_NOT_FINITE`); the checks run in this order. `qdt_insert_point*`
  rejects a point with a non-finite coordinate (`QDT_ERR_NOT_FINITE`) and one
  that would take the context beyond 100,000,000 points
  (`QDT_ERR_TOO_MANY_POINTS`); the context stays as it was. The geometric
  tests are exact for any overall scale of the input and for sets that mix very
  different magnitudes, with one exception: nonzero coordinates with an
  absolute value below about 2.2e-308 (subnormal numbers) have too few
  significant bits, and the output may then fail to be Delaunay without any
  error being reported.

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
  This takes time proportional to the size of the triangulation (about 0.69 s
  for 2.6 million points in the benchmark below, about 0.46 s with the parallel
  variant, which is about a third and about half of the time of the batch call
  itself) and the mesh stays in memory until the next batch call or
  `qdt_destroy`. The first constraint call after a batch call builds the same
  mesh.
- **Every later insertion** first locates the triangle that contains the new
  point, by walking through the triangulation from the triangle of the
  previous insertion, and then updates the triangulation around it. The cost
  therefore depends on how far the new point is from the previous one and
  grows only slowly with the total number of points: in the benchmark below,
  with uniformly random points, the median insertion takes about 0.005 ms at
  10,000 points, 0.027 ms at 100,000, 0.11 ms at 1 million and 0.18 ms at 2.6
  million points (40 times more for 260 times more points). A few insertions
  are much slower than the typical one: at 1 million and 2.6 million points the
  mean is 2.8 and 3.8 times the median. In the rare case that the walk cannot
  finish (it can only fail in a triangulation with constraints), the triangle
  is found by scanning all triangles, which takes time proportional to the
  number of triangles.
- In the benchmark below, a full rebuild with `qdt_execute` is faster than
  inserting points one by one once more than roughly 800 to 1,800 points
  (roughly 700 to 1,200 with the parallel variants, depending on the size)
  have to be added. These break-even figures use the mean insertion time, so
  they are on the cautious side; with the median they would be higher. If the
  point set changes by more than that, call `qdt_execute` or
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
  was, when `v0` or `v1` is not a point of the context
  (`QDT_ERR_INVALID_INDEX`), when `v0` and `v1` are the same point
  (`QDT_ERR_SAME_POINT`), when the context has no triangles
  (`QDT_ERR_NO_TRIANGLES`), when it would cross an edge that is already
  constrained (`QDT_ERR_CROSSES_CONSTRAINT`), or when the segment passes
  exactly through a third point (`QDT_ERR_THROUGH_POINT`; insert the pieces
  between consecutive points instead). The reasons are checked in this order;
  if a segment both crosses a constrained edge and passes through a point,
  either of the last two may be reported.
- The order of `v0` and `v1` does not matter. The index of a duplicate point
  acts as the lowest index with the same coordinates.
- Constraints are **kept** by `qdt_insert_point*`: a point inserted exactly on
  a constrained edge splits it, and both halves stay constrained. They are
  **discarded** by the next `qdt_execute` or `qdt_execute_parallel`.
- A single `qdt_insert_constraint` call finds its two end points in a table
  that is built together with the internal mesh and kept up to date by every
  insertion, so the lookup takes constant time and does not grow with the
  number of points. The cost of a constraint grows with the number of
  triangles its segment crosses. `qdt_insert_constraints_batch` gives the same
  result as a loop of single calls and reports the status of every entry; in
  the benchmark below it takes about the same time as the loop (within about
  2%).

```c
/* A flat quadrilateral: the Delaunay triangulation uses the diagonal 1-3.
   The constraint 0-2 forces the other diagonal. */
qdt_vec2d pts[] = { {0.0, 0.0}, {4.0, 0.0}, {3.0, 1.0}, {1.0, 1.0} };

qdt_t* q = qdt_create();
qdt_execute(q, pts, 4);          /* triangles: (3, 0, 1) and (1, 2, 3) */
qdt_insert_constraint(q, 0, 2);  /* returns QDT_OK; triangles: (2, 0, 1) and (0, 2, 3) */
qdt_insert_constraint(q, 1, 3);  /* returns QDT_ERR_CROSSES_CONSTRAINT; nothing changes */
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
  tiny inputs. In the benchmark below it takes about 4% longer than
  `qdt_execute` at 100,000 points, about 10% longer at 10,000 points, and much
  longer in relative terms for tiny inputs (about 15 times at 100 points). Use
  `qdt_execute` for such inputs.
- If the parallel implementation cannot run (its pre-allocated edge storage is
  exhausted, or the worker threads cannot be started), the call falls back to
  the serial algorithm and still returns the complete result with `QDT_OK`.

`qdt_insert_point_parallel` gives the same result as `qdt_insert_point`. The
insertion itself is single-threaded; only the one-time mesh build performed by
the first insertion after a batch call is parallel (and only when the
triangulation has at least 50,000 triangles). This matches the benchmark
below: `insert(0)` is faster in the parallel run, while later insertions take
about the same time in both.

### Error handling

Every function that can fail returns a `qdt_status` (see [Status
codes](#status-codes)); `qdt_status_string` turns a code into a text for log
messages. No function lets a C++ exception reach the caller: every function of
`qdt.h` catches all exceptions at the API boundary (an allocation failure
surfaces inside the library as `std::bad_alloc`). The functions can therefore
be called safely from C and from any other language that can call C functions.

What a failed call does to the context depends on the error:

- **The context is left exactly as it was** for `QDT_ERR_NULL_CONTEXT` (from
  any function), for `QDT_ERR_NOT_FINITE` and `QDT_ERR_TOO_MANY_POINTS` from
  `qdt_insert_point*`, for `QDT_ERR_NULL_ARGUMENT` from
  `qdt_insert_constraints_batch`, and for every rejection of a constraint
  (`QDT_ERR_INVALID_INDEX` up to `QDT_ERR_THROUGH_POINT`).
- **The context is emptied** (no points, no triangles, no constraints and no
  internal mesh) for `QDT_ERR_TOO_MANY_POINTS`, `QDT_ERR_NULL_ARGUMENT` and
  `QDT_ERR_NOT_FINITE` from `qdt_execute` and `qdt_execute_parallel`, which
  always replace the previous contents, and for `QDT_ERR_OUT_OF_MEMORY` and
  `QDT_ERR_INTERNAL` from any function. QDT discards the data instead of
  keeping it, because a failed operation may have updated some internal
  containers and not others, so what was left would not be a consistent
  triangulation. An emptied context stays valid: it can be queried, reused (for
  example with `qdt_execute`) and destroyed. `qdt_points_size()` returns `0`
  afterwards.
- `qdt_create` returns `NULL` if the memory for the context cannot be
  allocated.
- If `qdt_insert_constraints_batch` runs out of memory or hits an internal
  error, it stops, the context is emptied, and no entry counts as accepted;
  every entry of `results` then receives that code.
- The fallback of `qdt_execute_parallel` to the serial algorithm (exhausted
  edge storage, threads that cannot be started) is not an error: the call
  still returns the complete result with `QDT_OK`.

### Thread safety

A single context must not be used from several threads at the same time.
Different contexts share no state and can be used concurrently.

### Resource use

Measured with the benchmark below (Windows, committed private bytes per input
point, 100,000 to 2.6 million points):

- A context holds about 72 bytes per input point after a batch call: 16 for
  the point copy, 8 for the index map, and 48 for the triangles (two per
  point). Small inputs show somewhat more, for example 96 bytes per point at
  100,000 points.
- At its peak, a batch call uses about 150 bytes per point above the level
  before the call (about 200 with `qdt_execute_parallel`); about 72 of these
  remain afterwards.
- The mesh built by the first insertion or constraint call stays in memory and
  adds about 400 bytes per point (committed memory, including the table that
  maps point indices to mesh vertices). While it is being built, the usage
  peaks at about 520 bytes per point above the level after the batch call. The
  committed figure counts storage that is reserved but not fully used. On
  Linux, where only touched pages count, the resident set grew by about 180
  bytes per point for the mesh in a measurement with 100,000 points. Call
  `qdt_execute*` again or `qdt_destroy` to release the mesh.
- Internal indices are 32-bit, so inputs of more than 100,000,000 points are
  rejected (see [Results](#results)).

## Example

The complete, compilable program is [`example/example.c`](example/example.c).
These are the calls it makes, in order (`pts` is an array of `n` points; the
program also checks every returned `qdt_status`):

```c
qdt_t* q = qdt_create();
size_t index;

/* Batch-triangulate most of the points. */
qdt_execute(q, pts, n - 2);

/* Insert the remaining points incrementally. The last argument receives the
   index of the new point; it may be NULL if it is not needed. */
qdt_insert_point(q, pts[n - 2], &index);
qdt_insert_point(q, pts[n - 1], &index);

/* Use qdt_execute_parallel(q, pts, n) for a full parallel batch. */
/* Use qdt_insert_point_parallel(q, pt, &index) for parallel insertion. */

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
Point 4 inserted
Point 5 inserted
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

- Date: 2026-11-16 01:45 UTC
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
| 100       |       185 |      0.0323 |        0.4719 |   0.07x |            3.09 |              0.21 |
| 1,000     |     1,976 |      0.5788 |         1.012 |   0.57x |            1.73 |              0.99 |
| 10,000    |    19,968 |       6.237 |         6.871 |   0.91x |            1.60 |              1.46 |
| 100,000   |   199,958 |       66.89 |         69.74 |   0.96x |            1.49 |              1.43 |
| 1,000,000 | 1,999,963 |       742.6 |         508.2 |   1.46x |            1.35 |              1.97 |
| 2,624,409 | 5,248,777 |      1968.4 |         906.4 |   2.17x |            1.33 |              2.90 |

### Batch triangulation: point distributions

1,000,000 points per case; medians of 5 runs. Triangle counts differ between cases because the number of convex hull points and of duplicate points differs.

| Distribution             | Triangles | Serial (ms) | Parallel (ms) | Speedup | Serial (Mpts/s) |
| ------------------------ | --------: | ----------: | ------------: | ------: | --------------: |
| uniform random           | 1,999,953 |       742.0 |         526.6 |   1.41x |            1.35 |
| uniform, presorted by x  | 1,999,959 |       658.3 |         448.3 |   1.47x |            1.52 |
| gaussian                 | 1,999,980 |       759.8 |         550.9 |   1.38x |            1.32 |
| clustered (20 gaussians) | 1,999,977 |       760.0 |         535.5 |   1.42x |            1.32 |
| regular grid (shuffled)  | 1,996,002 |      1194.8 |         735.0 |   1.63x |            0.84 |
| jittered grid (shuffled) | 1,999,961 |       727.7 |         517.7 |   1.41x |            1.37 |
| all points on a circle   | 1,000,017 |       590.1 |         425.7 |   1.39x |            1.69 |
| flat band (aspect 1e6:1) | 1,999,962 |      1731.2 |        1130.7 |   1.53x |            0.58 |
| 50% duplicate points     |   999,964 |       434.5 |         453.3 |   0.96x |            2.30 |

### Incremental insertion

A base triangulation of N uniform random points, followed by 50 single point insertions. The first insertion builds the internal mesh and is much slower than the later ones. Break-even is the estimated number of inserted points at which a full rebuild becomes cheaper: 1 + (build - first insert) / mean(later inserts). Serial = `qdt_execute` + `qdt_insert_point`; parallel = `qdt_execute_parallel` + `qdt_insert_point_parallel`. Medians of up to 3 runs.

| Base points | Mode     | Build (ms) | 1st insert (ms) | Later median (ms) | Later mean (ms) | Later p95 (ms) | Break-even (points) |
| ----------- | -------- | ---------: | --------------: | ----------------: | --------------: | -------------: | ------------------: |
| 10,000      | serial   |      6.332 |           1.743 |            0.0046 |          0.0053 |         0.0090 |                ~871 |
| 10,000      | parallel |      6.920 |           1.767 |            0.0047 |          0.0056 |         0.0094 |                ~928 |
| 100,000     | serial   |      67.39 |           20.70 |            0.0265 |          0.0563 |         0.0727 |                ~830 |
| 100,000     | parallel |      70.52 |           14.88 |            0.0269 |          0.0579 |         0.0749 |                ~963 |
| 1,000,000   | serial   |      740.4 |           236.4 |            0.1093 |          0.3110 |         0.1855 |               ~1622 |
| 1,000,000   | parallel |      535.0 |           160.8 |            0.1109 |          0.3102 |         0.1895 |               ~1207 |
| 2,624,409   | serial   |     1982.1 |           686.7 |            0.1848 |          0.7082 |         0.2838 |               ~1830 |
| 2,624,409   | parallel |      948.3 |           460.4 |            0.1803 |          0.7285 |         0.2800 |                ~671 |

#### Per-insert detail, non-parallel, base 2,624,409 points

| Operation               | Triangles | Time (ms) |
| ----------------------- | --------: | --------: |
| Full (2,624,419 points) | 5,248,798 |    1983.8 |
| Most, 10 points left    | 5,248,778 |    1982.1 |
| insert(0)               | 5,248,780 |     689.9 |
| insert(1)               | 5,248,782 |    0.1459 |
| insert(2)               | 5,248,784 |    0.1988 |
| insert(3)               | 5,248,786 |    0.2416 |
| insert(4)               | 5,248,788 |    0.1999 |
| insert(5)               | 5,248,790 |    0.1023 |
| insert(6)               | 5,248,792 |    0.2131 |
| insert(7)               | 5,248,794 |    0.0804 |
| insert(8)               | 5,248,796 |    0.1548 |
| insert(9)               | 5,248,798 |    0.1891 |

#### Per-insert detail, parallel, base 2,624,409 points

| Operation               | Triangles | Time (ms) |
| ----------------------- | --------: | --------: |
| Full (2,624,419 points) | 5,248,798 |     939.3 |
| Most, 10 points left    | 5,248,778 |     931.8 |
| insert(0)               | 5,248,780 |     440.3 |
| insert(1)               | 5,248,782 |    0.1529 |
| insert(2)               | 5,248,784 |    0.2087 |
| insert(3)               | 5,248,786 |    0.2421 |
| insert(4)               | 5,248,788 |    0.2002 |
| insert(5)               | 5,248,790 |    0.0987 |
| insert(6)               | 5,248,792 |    0.2193 |
| insert(7)               | 5,248,794 |    0.0703 |
| insert(8)               | 5,248,796 |    0.1527 |
| insert(9)               | 5,248,798 |    0.1779 |

### Constrained Delaunay triangulation

A base triangulation of N uniform random points, then constraints from three sets: "short" (segments between nearby points, none of them crossing; some are already Delaunay edges), "long chords" (segments spanning the whole domain, each crossing about sqrt(N) triangles) and "ring" (a closed polygon through the points of a thin annulus). The first constraint is inserted with `qdt_insert_constraint` and includes the one-time internal mesh build; the others are inserted either with one `qdt_insert_constraints_batch` call or with a loop of `qdt_insert_constraint` calls. Every constraint of these sets is valid and must be accepted. Medians of up to 3 runs.

| Points    | Constraint set | Constraints | 1st call + mesh build (ms) | Batch (ms) | Batch (us/constraint) | Single calls (ms) | Single (us/constraint) |        Accepted |
| --------- | -------------- | ----------: | -------------------------: | ---------: | --------------------: | ----------------: | ---------------------: | --------------: |
| 10,000    | short          |         492 |                      1.537 |     0.5925 |                  1.21 |            0.5997 |                   1.22 |       492 / 492 |
| 10,000    | long chords    |          35 |                      1.662 |      2.417 |                  71.1 |             2.365 |                   69.6 |         35 / 35 |
| 10,000    | ring           |       1,842 |                      1.584 |      2.008 |                  1.09 |             1.978 |                   1.07 |   1,842 / 1,842 |
| 100,000   | short          |       4,839 |                      17.87 |      6.954 |                  1.44 |             6.966 |                   1.44 |   4,839 / 4,839 |
| 100,000   | long chords    |         100 |                      17.15 |      23.17 |                   234 |             23.22 |                    235 |       100 / 100 |
| 100,000   | ring           |       1,942 |                      16.90 |     0.9580 |                  0.49 |            0.9545 |                   0.49 |   1,942 / 1,942 |
| 1,000,000 | short          |      20,000 |                      203.7 |      33.94 |                  1.70 |             33.53 |                   1.68 | 20,000 / 20,000 |
| 1,000,000 | long chords    |         100 |                      201.5 |      94.67 |                   956 |             94.83 |                    958 |       100 / 100 |
| 1,000,000 | ring           |       2,006 |                      200.8 |     0.9529 |                  0.48 |            0.9644 |                   0.48 |   2,006 / 2,006 |

Accepted is shown as batch / single-call counts.

### Reading the result

`qdt_get_triangle` over all triangles of a uniform random point set.

| Points    | Triangles | Time (ms) | ns per call | Mcalls/s |
| --------- | --------: | --------: | ----------: | -------: |
| 1,000,000 | 1,999,957 |     5.275 |         2.6 |    379.2 |

### Independent contexts on several threads

Every thread triangulates the same 100,000 uniform random points with `qdt_execute` on its own context (the contexts share no state). Scaling is the aggregate throughput divided by the throughput of one thread; efficiency is scaling divided by the thread count. Hardware threads: 12. Medians of 5 runs.

| Threads | Wall time (ms) | Aggregate (Mpts/s) | Scaling | Efficiency |
| ------- | -------------: | -----------------: | ------: | ---------: |
| 1       |          68.02 |               1.47 |   1.00x |       100% |
| 2       |          69.72 |               2.87 |   1.95x |        98% |
| 4       |          80.19 |               4.99 |   3.39x |        85% |
| 8       |          103.4 |               7.74 |   5.26x |        66% |
| 12      |          118.0 |              10.17 |   6.92x |        58% |

### Memory use

Committed private bytes (PrivateUsage) measured around the library calls, in a freshly started process per measurement; the input point array is part of the baseline. "Retained" is the growth after the call. "Peak" is the highest level during the call minus the level before it. Mesh = the internal half-edge mesh built by the first insertion after a batch call. Bytes are per input point.

| Points    | Mode     | Retained after batch (B/pt) | Peak during batch (B/pt) | Mesh retained (B/pt) | Peak during mesh build (B/pt) | Retained, batch + mesh (MiB) |
| --------- | -------- | --------------------------: | -----------------------: | -------------------: | ----------------------------: | ---------------------------: |
| 100,000   | serial   |                        96.4 |                    152.5 |                401.4 |                         509.7 |                         47.5 |
| 100,000   | parallel |                        96.3 |                    216.4 |                409.8 |                         518.1 |                         48.3 |
| 1,000,000 | serial   |                        72.2 |                    152.3 |                402.8 |                         523.1 |                        453.0 |
| 1,000,000 | parallel |                        72.2 |                    202.5 |                402.8 |                         527.1 |                        453.1 |
| 2,624,409 | serial   |                        72.1 |                    152.3 |                402.8 |                         519.0 |                       1188.7 |
| 2,624,409 | parallel |                        72.2 |                    201.2 |                402.8 |                         519.0 |                       1188.7 |

### Comparison with CGAL

CGAL 6.2.1 `Delaunay_triangulation_2` with `Exact_predicates_inexact_constructions_kernel`, inserting the whole point range at once (CGAL sorts spatially itself). Point conversion is not timed. Uniform random points, medians of 5 runs (3 for the largest sizes). The speedup columns are CGAL time divided by QDT time. Make sure CGAL was compiled with NDEBUG.

| Points    | CGAL (ms) | QDT serial (ms) | QDT parallel (ms) | Serial speedup | Parallel speedup | Triangles (CGAL) | Triangles (QDT) |
| --------- | --------: | --------------: | ----------------: | -------------: | ---------------: | ---------------: | --------------: |
| 10,000    |     6.186 |           6.157 |             6.778 |          1.00x |            0.91x |           19,968 |          19,968 |
| 100,000   |     65.91 |           66.59 |             69.93 |          0.99x |            0.94x |          199,958 |         199,958 |
| 1,000,000 |     696.4 |           738.5 |             514.1 |          0.94x |            1.35x |        1,999,963 |       1,999,963 |
| 2,624,409 |    1907.5 |          1972.1 |             948.8 |          0.97x |            2.01x |        5,248,777 |       5,248,777 |

### Notes on the results

- **Batch.** Serial throughput falls from about 3.1 million points per second
  for tiny inputs to about 1.3 to 1.4 million for millions of points.
  2,624,409 points (5,248,777 triangles) take 1.97 s with `qdt_execute` and
  0.91 s with `qdt_execute_parallel` (2.17x). The parallel call only pays off
  for large inputs; see [Parallel mode](#parallel-mode).
- **Distributions.** Compared with uniform random points (742 ms for 1,000,000
  points), presorted input is 11% faster and points on a circle are 20%
  faster; gaussian, clustered, and jittered-grid inputs are within 3%; a
  regular grid with many cocircular points is 1.6x slower, and the very thin
  band is 2.3x slower. The input with 50% duplicates is the fastest in absolute
  terms because only 500,000 distinct points are triangulated. That is below
  the 524,288-point threshold for splitting work between threads, so it gets
  no parallel speedup (0.96x).
- **Incremental insertion.** The first insertion costs about a third of a
  serial build (687 ms against 1,982 ms at 2,624,409 points) and about half of
  a parallel build (460 ms against 948 ms). Later insertions grow slowly with
  the size: the median is 0.005 ms at 10,000 points, 0.027 ms at 100,000, 0.11
  ms at 1,000,000 and 0.18 ms at 2,624,409. At the larger sizes a few
  insertions are much slower than the typical one, so the mean (0.31 ms and
  0.71 ms at 1,000,000 and 2,624,409 points) is well above the median; the
  break-even figures are computed from the mean. Building most of the points
  and then inserting the last 10 one by one takes 1.35x (serial) or 1.46x
  (parallel) as long as triangulating all of them at once, because the first
  insertion alone costs a third or half of a build. Incremental insertion
  therefore pays off only for a small number of points: a full rebuild is
  cheaper once roughly 700 to 1,800 points are added.
- **Constraints.** At 1,000,000 points the first constraint call, which
  includes building the internal mesh, takes about 200 ms (27% of the 743 ms
  batch triangulation). After that, constraints between nearby points cost
  about 1.7 microseconds each at 1,000,000 points and about 0.5 microseconds
  each for the thin ring. A single `qdt_insert_constraint` call takes as long as
  the same constraint inside `qdt_insert_constraints_batch` (within about 2% in
  every row above), because the end points are found through a lookup table.
  Long chords, each crossing about sqrt(N) triangles, cost about 0.07 ms at
  10,000 points, 0.23 ms at 100,000 and 0.96 ms at 1,000,000 points; the
  crossing dominates there.
- **Several contexts.** Independent contexts scale well: 3.4x on 4 threads and
  6.9x on all 12 hardware threads of this 6-core CPU (the second thread of a
  core adds less than a full core).
- **Reading the result.** `qdt_get_triangle` takes about 2.6 ns per call.
- **CGAL.** Both libraries produce the same number of triangles at every size.
  Serial QDT and CGAL take about the same time: QDT needs 1.00x, 1.01x, 1.06x
  and 1.03x of CGAL's time at 10,000, 100,000, 1,000,000 and 2,624,409 points,
  so CGAL is up to 6% faster in serial mode. CGAL's 2D triangulation runs on
  the calling thread. `qdt_execute_parallel` is slightly slower than CGAL while
  it still runs on one thread (0.91x at 10,000 and 0.94x at 100,000 points; see
  [Parallel mode](#parallel-mode)), and faster once it splits the work: 1.35x
  at 1,000,000 points and 2.01x at 2,624,409 points. The comparison covers only
  batch triangulation of uniform random points on one machine and compiler;
  CGAL also offers other kernels, other point orders and other operations
  that this benchmark does not measure.

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
small sizes only, and `--scenario LIST` selects scenarios (`batch`, `dist`,
`insert`, `constraint`, `query`, `concurrent`, `memory`, `cgal`). The
first-insertion time varies the most between runs.

The `cgal` scenario is optional. It is compiled in only with
`-DQDT_BENCH_CGAL`, needs CGAL 5.x or newer (and therefore C++17), and is
skipped silently when the macro is not defined. Build with NDEBUG so that CGAL
is not slowed down by its assertions. On Linux, with the CGAL headers and
GMP/MPFR installed (link GMP/MPFR only if your CGAL installation uses them):

```sh
g++ -std=c++17 -O2 -DNDEBUG -DQDT_BENCH_CGAL -pthread -Ilibqdt \
    test/qdt_bench.cpp libqdt.a -lgmp -lmpfr -lm -o qdt_bench
./qdt_bench --scenario cgal --output cgal.md
```

On Windows, use `/std:c++17` and `/DQDT_BENCH_CGAL` in the MSVC command above
and add the include and library paths of your CGAL installation.

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
  existing triangulation. The first insertion builds an internal mesh (about a
  third of the time of a serial batch build), and in the benchmark a full
  rebuild is cheaper once roughly 700 to 1,800 points have to be added. If the
  point set changes substantially, call `qdt_execute` or
  `qdt_execute_parallel` again with the complete updated point array.
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
- Errors are reported only as `qdt_status` codes. A function that runs out of
  memory or hits an internal error empties the context, and so does a rejected
  batch call (`qdt_execute*`); the rejections of single points and constraints
  leave it unchanged (see [Error handling](#error-handling)).
- Inputs must contain only finite coordinates, and a context must not be shared
  between threads without external synchronization.

## License

QDT is released under the MIT License. See [LICENSE](LICENSE).