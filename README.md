# QDT - Quick 2D Delaunay Triangulation

QDT is a small library for fast Delaunay triangulation of 2D point sets. It is
implemented in C++11 and exposes a plain C API built around an opaque handle,
so it can be used from C, C++, and any language that can call C functions. It
supports single-threaded and parallel batch triangulation as well as
incremental point insertion.

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
- Suitable for point clouds, mesh generation, interpolation, and related
  geometry processing tasks.
- MIT licensed.

## Requirements

Using the library:

- Any C99 (or later) or C++ compiler for code that includes `qdt.h`. The header
  only depends on `<stddef.h>`.
- On Windows, see [Windows notes](#windows-notes) for the `QDT_DLL` macro.

Building the library from source:

- A C++11 (or later) compiler and a threads library (`-pthread` on Linux). The
  sources compile with GCC 13 using `-std=c++11`.
- Source files: `qdt.cpp`, `GuibasStolfiDwyer_Wrapper.cpp`,
  `GuibasStolfiDwyer_Threads.cpp`, plus the headers `qdt.h` and `wrapper.hxx`.

## Building

Linux and macOS, shared library:

```sh
g++ -std=c++11 -O2 -pthread -fPIC -shared \
    qdt.cpp GuibasStolfiDwyer_Wrapper.cpp GuibasStolfiDwyer_Threads.cpp \
    -o libqdt.so
```

Static library:

```sh
g++ -std=c++11 -O2 -pthread -fPIC -c \
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

typedef struct qdt_t qdt_t; /* Opaque triangulation context. */
```

### Functions

| Function | Description |
| --- | --- |
| `qdt_t* qdt_create(void)` | Creates a new, empty triangulation context. See [Error handling](#error-handling). |
| `void qdt_execute(qdt_t* ptr, qdt_vec2d* pts, size_t n)` | Computes the Delaunay triangulation of `n` points in one batch. Replaces everything previously stored in the context, including inserted points. The points are copied; `pts` is not modified and need not stay valid after the call. |
| `void qdt_execute_parallel(qdt_t* ptr, qdt_vec2d* pts, size_t n)` | Same as `qdt_execute`, but may use multiple threads. See [Parallel mode](#parallel-mode). |
| `size_t qdt_insert_point(qdt_t* ptr, qdt_vec2d pt)` | Appends one point to the context and updates the triangulation incrementally. Returns the total number of points in the context after the call, so the new point's index is the return value minus 1. Returns `0` if `ptr` is `NULL`. |
| `size_t qdt_insert_point_parallel(qdt_t* ptr, qdt_vec2d pt)` | Same result as `qdt_insert_point`. Only a one-time setup step uses multiple threads. See [Parallel mode](#parallel-mode). |
| `size_t qdt_triangles_size(qdt_t* ptr)` | Returns the current number of triangles. |
| `size_t qdt_get_triangle(qdt_t* ptr, size_t T[3], size_t id)` | Writes the three point indices of triangle `id` to `T`. Returns `1` when the triangle exists, `0` otherwise (`T` is left untouched). |
| `void qdt_destroy(qdt_t* ptr)` | Destroys the context and releases its resources. `NULL` is allowed. |

Passing a `NULL` context to any function is safe: the functions do nothing and
return `0`.

### Results

- **Point indices.** For `qdt_execute` and `qdt_execute_parallel`, indices
  refer to the input array order. Points added with `qdt_insert_point` or
  `qdt_insert_point_parallel` are appended in insertion order. Indices stay
  valid until the next batch call.
- **Triangle ids.** Valid ids are `0 .. qdt_triangles_size() - 1`. Ids are
  **not stable**: any insertion can remove triangles and reuse their ids, so
  do not cache ids across `qdt_execute*` or `qdt_insert_point*` calls.
- **Orientation.** The vertices of each triangle are listed counter-clockwise
  (positive signed area in a right-handed coordinate system). Which vertex
  comes first is unspecified, as is the order of the triangles.
- **Duplicate points.** If several points have identical coordinates, only the
  one with the lowest index appears in the triangles. The others still occupy
  an index (and count towards the value returned by `qdt_insert_point`), but
  are never referenced by a triangle.
- **Degenerate input.** With fewer than three distinct points, or with all
  points collinear, the result has zero triangles.
- **Input requirements.** Coordinates must be finite. `NaN` and infinities are
  not supported.

### Batch and incremental usage

QDT supports two main workflows:

1. **Batch triangulation:** call `qdt_execute` or `qdt_execute_parallel` with
   the full point array.
2. **Incremental insertion:** build the triangulation for most points with
   `qdt_execute` or `qdt_execute_parallel`, then add the remaining points one
   by one with `qdt_insert_point` or `qdt_insert_point_parallel`.

Incremental insertion is meant for adding a small number of points to an
existing triangulation. Keep these costs in mind:

- The **first insertion** after a batch call builds an internal half-edge mesh.
  This takes time proportional to the size of the triangulation (about 0.4-0.5
  s for 2.6 million points in the benchmark below) and the mesh stays in memory
  until the next batch call or `qdt_destroy`.
- **Every insertion** costs time proportional to the current number of points,
  not to the size of the affected neighborhood. It is not a constant-time or
  logarithmic operation.
- In the benchmark below, a full rebuild with `qdt_execute` is faster than
  inserting points one by one once more than roughly 100 points (about 50 with
  the parallel variants) have to be added. If the point set changes by more
  than that, call `qdt_execute` or `qdt_execute_parallel` again with the
  complete point array.
- Inserting into a context that has no triangles yet (fewer than three
  distinct points, or all points collinear) works, but each insertion then
  triggers a full rebuild.

### Parallel mode

`qdt_execute_parallel` produces the same set of triangles as `qdt_execute`,
but not necessarily in the same order.

- The thread count is `std::thread::hardware_concurrency()`, capped at 64, and
  is not configurable. Threads exist only for the duration of the call.
- With fewer than 200,000 distinct points (or a single hardware thread), the
  call runs on the calling thread. For inputs up to a few hundred thousand
  points, the parallel version may be no faster than `qdt_execute`.

`qdt_insert_point_parallel` gives the same result as `qdt_insert_point`. The
insertion itself is single-threaded; only the one-time mesh build performed by
the first insertion after a batch call is parallel (and only when the
triangulation has at least 50,000 triangles). This matches the benchmark
below: `insert(0)` is faster in the parallel run, while later insertions take
the same time in both.

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

Estimated from data-structure sizes, not measured:

- A context holds roughly 70 bytes per input point after a batch call (point
  copy, index map, triangles).
- A batch call needs roughly 100-150 additional bytes per point while it runs.
- The mesh built by the first insertion adds roughly 140 bytes per point.
- Internal indices are 32-bit, so inputs of around 100 million points or more
  are not supported.

## Example

```c
#include "qdt.h"
#include <stdio.h>

int main(void)
{
    /* Points 0-3 form a square (all four are cocircular). Point 4 is its
       center, and point 5 lies exactly on the edge between points 4 and 2. */
    qdt_vec2d pts[] = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},
        {1.0, 1.0},
        {0.5, 0.5},
        {0.25, 0.75}
    };
    size_t n = sizeof(pts) / sizeof(pts[0]);

    qdt_t* q = qdt_create();
    if (!q)
    {
        return 1;
    }

    /* Batch-triangulate most of the points. */
    qdt_execute(q, pts, n - 2);

    /* Insert the remaining points incrementally. */
    qdt_insert_point(q, pts[n - 2]);
    qdt_insert_point(q, pts[n - 1]);

    /* Use qdt_execute_parallel(q, pts, n) for a full parallel batch. */
    /* Use qdt_insert_point_parallel(q, pt) for parallel insertion. */

    size_t triangle_count = qdt_triangles_size(q);
    printf("Triangles: %zu\n", triangle_count);

    for (size_t i = 0; i < triangle_count; ++i)
    {
        size_t tri[3];
        if (qdt_get_triangle(q, tri, i))
        {
            printf("Triangle %zu: %zu, %zu, %zu\n",
                   i, tri[0], tri[1], tri[2]);
        }
    }

    qdt_destroy(q);
    return 0;
}
```

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

The following benchmark uses 2,624,409 input points. Point generation
(9685.83 ms) is not included in the triangulation timings below. Absolute
timings depend on CPU, compiler, optimization flags, and thread count. The
relative comparison is more useful.

### Batch comparison

```text
CGAL generates 2624409 points, 5248777 triangles and 7873185 edges in 2039.27 ms
QDT generates (Full) 5248777 triangles in 1833.38 ms
QDT (Parallel) generates (Full) 5248777 triangles in 951.955 ms
```

Relative performance, using CGAL as the baseline (`CGAL = 1.00`):

| Implementation | Time | Relative time (CGAL = 1.00) | Relative speed |
| --- | ---: | ---: | ---: |
| CGAL | 2039.27 ms | 1.00 | 1.00x |
| QDT (Full) | 1833.38 ms | 0.90 | 1.11x |
| QDT (Parallel, Full) | 951.955 ms | 0.47 | 2.14x |

In this test, single-threaded QDT takes about 0.90x the time of CGAL, or
roughly 1.11x the speed. Parallel QDT takes about 0.47x the time of CGAL, or
roughly 2.14x the speed. Parallel QDT is about 1.93x faster than
single-threaded QDT. These ratios are approximate and may vary with hardware
and workload.

### Incremental insertion

In this test, `Full` means all points are triangulated in one batch. `Most`
means most points are triangulated first, leaving 10 points to be inserted
one by one. `insert(i)` means inserting the `i`-th remaining point. After all
10 insertions, the triangle count reaches the same value as `Full`.

#### Non-parallel

| Operation | Triangles | Time |
| --- | ---: | ---: |
| Full | 5,248,777 | 1833.38 ms |
| Most, 10 points left | 5,248,758 | 1831.38 ms |
| insert(0) | 5,248,760 | 529.031 ms |
| insert(1) | 5,248,762 | 11.5622 ms |
| insert(2) | 5,248,764 | 11.4107 ms |
| insert(3) | 5,248,766 | 11.4971 ms |
| insert(4) | 5,248,768 | 12.311 ms |
| insert(5) | 5,248,769 | 11.6095 ms |
| insert(6) | 5,248,772 | 11.2689 ms |
| insert(7) | 5,248,773 | 11.8612 ms |
| insert(8) | 5,248,774 | 11.3663 ms |
| insert(9) | 5,248,777 | 18.0396 ms |

#### Parallel

| Operation | Triangles | Time |
| --- | ---: | ---: |
| Full | 5,248,777 | 951.955 ms |
| Most, 10 points left | 5,248,758 | 917.471 ms |
| insert(0) | 5,248,760 | 375.75 ms |
| insert(1) | 5,248,762 | 11.0857 ms |
| insert(2) | 5,248,764 | 10.7463 ms |
| insert(3) | 5,248,766 | 10.6544 ms |
| insert(4) | 5,248,768 | 11.9406 ms |
| insert(5) | 5,248,769 | 10.8896 ms |
| insert(6) | 5,248,772 | 10.7014 ms |
| insert(7) | 5,248,773 | 10.8919 ms |
| insert(8) | 5,248,774 | 10.6272 ms |
| insert(9) | 5,248,777 | 10.7781 ms |

The first insertion is noticeably slower than later ones because it builds the
internal mesh (parallelized in the parallel run). After that, each insertion
takes roughly 10-18 ms in the non-parallel test and roughly 10.6-11.9 ms in the
parallel test; the two are the same code path, so the difference is noise.
Per-insertion time grows with the number of points in the context.

`Most + 10 insertions` takes more total time than `Full` in this benchmark, so
incremental insertion is not a replacement for full batch triangulation when
all points are known up front. Its advantage is avoiding a full rebuild when
only a small number of new points need to be added to an existing
triangulation. At this problem size, a rebuild (1833 ms single-threaded, 952 ms
parallel) wins once roughly 100 (single-threaded) or 50 (parallel) points have
to be added.

## Current limitations

- Incremental insertion is intended for adding a small number of points to an
  existing triangulation. Its per-point cost is proportional to the current
  number of points. If the point set changes substantially, call
  `qdt_execute` or `qdt_execute_parallel` again with the complete updated
  point array.
- The current API does not provide deletion of points or direct editing of
  existing vertices.
- Triangle ids are not stable across insertions, and the order of triangles
  and of the vertices within a triangle is unspecified. Point indices are
  stable; keep the input point order and record insertion order if you need
  to map results back to your data.
- Parallel functions may use multiple threads, but actual speedup depends on
  hardware, compiler, thread count, and workload size. The thread count cannot
  be configured.
- Errors are reported as C++ exceptions, not return codes (see
  [Error handling](#error-handling)).
- Inputs must contain only finite coordinates, and a context must not be shared
  between threads without external synchronization.

## License

QDT is released under the MIT License. See [LICENSE](LICENSE).