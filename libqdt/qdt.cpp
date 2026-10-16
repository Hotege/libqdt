
#include "qdt.h"
#include "wrapper.hxx"

namespace wrapper
{
    using Point = qdt_vec2d;
    Triangulation delaunay_triangulate(Point* pts, size_t sz);
}

namespace parallel
{
    using Point = qdt_vec2d;
    wrapper::Triangulation delaunay_triangulate(Point* pts, size_t sz);
}

qdt_t* qdt_create()
{
    return reinterpret_cast<qdt_t*>(new wrapper::Triangulation);
}

void qdt_execute(qdt_t* ptr, qdt_vec2d* pts, size_t n)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return;
    *Q = wrapper::delaunay_triangulate(pts, n);
}

void qdt_execute_parallel(qdt_t* ptr, qdt_vec2d* pts, size_t n)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return;
    *Q = parallel::delaunay_triangulate(pts, n);
}

size_t qdt_insert_point(qdt_t* ptr, const qdt_vec2d& pt)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return 0;
    Q->insert(pt);
    return Q->input_points.size();
}

size_t qdt_insert_point_parallel(qdt_t* ptr, const qdt_vec2d& pt)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return 0;
    Q->insert_parallel(pt);
    return Q->input_points.size();
}

size_t qdt_triangles_size(qdt_t* ptr)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (Q)
        return Q->triangles.size();
    return 0;
}

size_t qdt_get_triangle(qdt_t* ptr, size_t T[3], size_t id)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return 0;
    if (id >= Q->triangles.size())
        return 0;
    for (int i = 0; i < 3; i++)
        T[i] = Q->triangles[id][i];
    return 1;
}

void qdt_destroy(qdt_t* ptr)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    delete Q;
}
