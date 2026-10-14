# QDT - Quick 2D Delaunay Triangulation

QDT is a small C library for quick Delaunay triangulation of 2D point sets.
It exposes a simple C API and supports both single-threaded and parallel
execution.

## Features

- Quick 2D Delaunay triangulation.
- Simple C API with an opaque handle.
- Single-threaded and parallel execution modes.
- Batch triangulation for a complete point set. Incremental insertion is not
  currently supported.
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
| `void qdt_execute(qdt_t* ptr, qdt_vec2d* pts, size_t n)` | Computes the Delaunay triangulation for `n` points. |
| `void qdt_execute_parallel(qdt_t* ptr, qdt_vec2d* pts, size_t n)` | Same as `qdt_execute`, but may use multiple threads. |
| `size_t qdt_triangles_size(qdt_t* ptr)` | Returns the number of triangles from the last execution. |
| `void qdt_get_triangle(qdt_t* ptr, size_t T[3], size_t id)` | Gets the three point indices of triangle `id`. |
| `void qdt_destroy(qdt_t* ptr)` | Destroys the context and releases resources. |

Triangle vertex indices returned by `qdt_get_triangle` refer to the input
point array passed to `qdt_execute` or `qdt_execute_parallel`.

**Batch-only API:** `qdt_execute` and `qdt_execute_parallel` compute the
triangulation for the entire input point set in one call. QDT currently does
not provide an incremental insertion API for adding points to an existing
triangulation. To update the triangulation for a modified point set, call
`qdt_execute` or `qdt_execute_parallel` again with the full point array.

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
        {1.0, 1.0}
    };
    size_t n = sizeof(pts) / sizeof(pts[0]);

    qdt_t* q = qdt_create();
    if (!q)
    {
        return 1;
    }

    /* Use qdt_execute_parallel(q, pts, n); for large point sets. */
    qdt_execute(q, pts, n);

    size_t triangle_count = qdt_triangles_size(q);
    printf("Triangles: %zu\n", triangle_count);

    for (size_t i = 0; i < triangle_count; ++i)
    {
        size_t tri[3];
        qdt_get_triangle(q, tri, i);
        printf("Triangle %zu: %zu, %zu, %zu\n",
               i, tri[0], tri[1], tri[2]);
    }

    qdt_destroy(q);
    return 0;
}
```

## Performance

The following benchmark was run on a local machine with 2,624,409 input
points. Absolute timings depend on CPU, compiler, optimization flags, and
thread count. The relative comparison is more useful.

```text
CGAL generates 2624409 points, 5248777 triangles and 7873185 edges in 2039.27 ms
QDT generates 5248777 triangles in 1813.65 ms
QDT (Parallel) generates 5248777 triangles in 899.179 ms
```

Relative performance, using CGAL as the baseline (`CGAL = 1.00`):

| Implementation | Time | Relative time (CGAL = 1.00) | Relative speed |
| --- | ---: | ---: | ---: |
| CGAL | 2039.27 ms | 1.00 | 1.00x |
| QDT | 1813.65 ms | 0.89 | 1.12x |
| QDT (Parallel) | 899.179 ms | 0.44 | 2.27x |

In this test, single-threaded QDT takes about 0.89x the time of CGAL, or
roughly 1.12x the speed. Parallel QDT takes about 0.44x the time of CGAL,
or roughly 2.27x the speed. Parallel QDT is about 2.02x faster than
single-threaded QDT. These ratios are approximate and may vary with hardware
and workload.

## Current limitations

- QDT currently supports one-shot batch triangulation only. After
  `qdt_execute` or `qdt_execute_parallel`, the triangulation is built for the
  entire input array.
- There is no incremental insertion API for adding new points to an existing
  triangulation. If the point set changes, call `qdt_execute` or
  `qdt_execute_parallel` again with the complete updated point array.

## Building

Add the directory containing `qdt.h` to your include path and link against
the QDT library, for example `libqdt.a`, `libqdt.so`, or `qdt.dll`.

On Windows, define `QDT_DLL` when building the DLL:

```c
#define QDT_DLL
```