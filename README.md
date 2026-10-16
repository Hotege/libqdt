# QDT - Quick 2D Delaunay Triangulation

QDT is a small C library for quick Delaunay triangulation of 2D point sets.
It exposes a simple C API and supports both single-threaded and parallel
execution, as well as batch triangulation and incremental point insertion.

## Features

- Quick 2D Delaunay triangulation.
- Simple C API with an opaque handle.
- Single-threaded and parallel execution modes.
- Batch triangulation for a complete point set.
- Incremental insertion of individual points, with single-threaded and
  parallel variants.
- Suitable for point clouds, mesh generation, interpolation, and related
  geometry processing tasks.
- MIT licensed.

## Requirements

- A C compiler with `size_t` support, such as C99 or later.
- The QDT library or source files linked into your project.
- On Windows, define `QDT_DLL` when building the DLL. Consumers that link
  against the DLL use the default `dllimport` behavior.

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
| `qdt_t* qdt_create()` | Creates a new triangulation context. Returns `NULL` on allocation failure. |
| `void qdt_execute(qdt_t* ptr, qdt_vec2d* pts, size_t n)` | Computes the Delaunay triangulation for `n` points in one batch. |
| `void qdt_execute_parallel(qdt_t* ptr, qdt_vec2d* pts, size_t n)` | Same as `qdt_execute`, but may use multiple threads. |
| `size_t qdt_insert_point(qdt_t* ptr, qdt_vec2d pt)` | Inserts one point into the current triangulation. Returns the total number of points after insertion. |
| `size_t qdt_insert_point_parallel(qdt_t* ptr, qdt_vec2d pt)` | Same as `qdt_insert_point`, but may use multiple threads. |
| `size_t qdt_triangles_size(qdt_t* ptr)` | Returns the number of triangles from the last execution or insertion. |
| `size_t qdt_get_triangle(qdt_t* ptr, size_t T[3], size_t id)` | Gets the three point indices of triangle `id`. Returns `1` when found, `0` otherwise. |
| `void qdt_destroy(qdt_t* ptr)` | Destroys the context and releases resources. |

Triangle vertex indices returned by `qdt_get_triangle` refer to points in the
current triangulation context. For `qdt_execute` and `qdt_execute_parallel`,
they correspond to the input point array order. For `qdt_insert_point` and
`qdt_insert_point_parallel`, inserted points are appended in insertion order;
the index of a newly inserted point is `qdt_insert_point(...) - 1`.

### Batch and incremental usage

QDT supports two main workflows:

1. **Batch triangulation:** call `qdt_execute` or `qdt_execute_parallel` with
   the full point array.
2. **Incremental insertion:** build the triangulation for most points with
   `qdt_execute` or `qdt_execute_parallel`, then add remaining points one by
   one with `qdt_insert_point` or `qdt_insert_point_parallel`.

Incremental insertion is useful when only a small number of points need to be
added to an existing triangulation. If the point set changes substantially,
calling `qdt_execute` or `qdt_execute_parallel` again with the complete point
array may be simpler and faster.

## Example

```c
#include "qdt.h"
#include <stdio.h>

int main(void)
{
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

## Performance

The following benchmark uses 2,624,409 input points. In this run, generating
the points took 9685.83 ms; point generation time is not included in the
triangulation timings below. Absolute timings depend on CPU, compiler,
optimization flags, and thread count. The relative comparison is more useful.

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

The first insertion is noticeably slower than later ones. After that, each
insertion takes roughly 10-18 ms in the non-parallel test and roughly
10.6-11.9 ms in the parallel test. Note that `Most + 10 insertions` takes
more total time than `Full` in this benchmark, so incremental insertion is
not a replacement for full batch triangulation when all points are known up
front. Its advantage is avoiding a full rebuild when only a small number of
new points need to be added to an existing triangulation.

## Current limitations

- Incremental insertion is intended for adding a small number of points to an
  existing triangulation. If the point set changes substantially, call
  `qdt_execute` or `qdt_execute_parallel` again with the complete updated
  point array.
- The current API does not provide deletion of points or direct editing of
  existing vertices.
- Parallel functions may use multiple threads, but actual speedup depends on
  hardware, compiler, thread count, and workload size.
- Triangle vertex order and indices may depend on the implementation. Keep the
  input point order and record insertion order if stable indices are needed.

## Building

Add the directory containing `qdt.h` to your include path and link against
the QDT library, for example `libqdt.a`, `libqdt.so`, or `qdt.dll`.

On Windows, define `QDT_DLL` when building the DLL:

```c
#define QDT_DLL
```