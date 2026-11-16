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
    qdt_status st;
    size_t index;

    qdt_t* q = qdt_create();
    if (!q)
    {
        return 1;
    }

    /* Batch-triangulate most of the points. */
    st = qdt_execute(q, pts, n - 2);
    if (st != QDT_OK)
    {
        fprintf(stderr, "qdt_execute: %s\n", qdt_status_string(st));
        qdt_destroy(q);
        return 1;
    }

    /* Insert the remaining points incrementally. The last argument receives
       the index of the new point; it may be NULL if it is not needed. */
    st = qdt_insert_point(q, pts[n - 2], &index);
    if (st != QDT_OK)
    {
        fprintf(stderr, "qdt_insert_point: %s\n", qdt_status_string(st));
        qdt_destroy(q);
        return 1;
    }
    printf("Point %zu inserted\n", index);

    st = qdt_insert_point(q, pts[n - 1], &index);
    if (st != QDT_OK)
    {
        fprintf(stderr, "qdt_insert_point: %s\n", qdt_status_string(st));
        qdt_destroy(q);
        return 1;
    }
    printf("Point %zu inserted\n", index);

    /* Use qdt_execute_parallel(q, pts, n) for a full parallel batch. */
    /* Use qdt_insert_point_parallel(q, pt, &index) for parallel insertion. */

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