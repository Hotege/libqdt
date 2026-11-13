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