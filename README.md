# QDT - Quick 2D Delaunay Triangulation

QDT is a small C library for quick Delaunay triangulation of 2D point sets.
It exposes a simple C API and supports both single-threaded and parallel
execution.

## Features

- Quick 2D Delaunay triangulation.
- Simple C API with an opaque handle.
- Single-threaded and parallel execution modes.
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

## Building

Add the directory containing `qdt.h` to your include path and link against
the QDT library, for example `libqdt.a`, `libqdt.so`, or `qdt.dll`.

On Windows, define `QDT_DLL` when building the DLL:

```c
#define QDT_DLL
```

## License

MIT License

Copyright (c) 2026 QDT contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
