#include <cmath>

#include "qdt.h"
#include "wrapper.hxx"

namespace wrapper
{
    using Point = qdt_vec2d;
    Triangulation delaunay_triangulate(const Point* pts, size_t sz);
}

namespace parallel
{
    using Point = qdt_vec2d;
    wrapper::Triangulation delaunay_triangulate(const Point* pts, size_t sz);
}

namespace
{
    bool point_ok(const qdt_vec2d& p)
    {
        return std::isfinite(p.x) && std::isfinite(p.y);
    }

    // Input that a batch call can process: at most wrapper::kMaxPoints points,
    // a non-NULL array unless it is empty, and finite coordinates only.
    bool batch_input_ok(const qdt_vec2d* pts, size_t n)
    {
        if (n > wrapper::kMaxPoints)
            return false;
        if (n == 0)
            return true;
        if (!pts)
            return false;
        for (size_t i = 0; i < n; ++i)
        {
            if (!point_ok(pts[i]))
                return false;
        }
        return true;
    }
}

qdt_t* qdt_create()
{
    return reinterpret_cast<qdt_t*>(new wrapper::Triangulation);
}

void qdt_execute(qdt_t* ptr, const qdt_vec2d* pts, size_t n)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return;
    if (!batch_input_ok(pts, n))
    {
        // Rejected input: the call still replaces the previous contents, so
        // the context is left empty rather than holding stale results.
        *Q = wrapper::Triangulation();
        return;
    }
    *Q = wrapper::delaunay_triangulate(pts, n);
}

void qdt_execute_parallel(qdt_t* ptr, const qdt_vec2d* pts, size_t n)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return;
    if (!batch_input_ok(pts, n))
    {
        *Q = wrapper::Triangulation();
        return;
    }
    *Q = parallel::delaunay_triangulate(pts, n);
}

size_t qdt_insert_point(qdt_t* ptr, qdt_vec2d pt)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return 0;
    // A rejected point is not added: the unchanged point count is returned.
    if (!point_ok(pt) || Q->input_points.size() >= wrapper::kMaxPoints)
        return Q->input_points.size();
    Q->insert(pt);
    return Q->input_points.size();
}

size_t qdt_insert_point_parallel(qdt_t* ptr, qdt_vec2d pt)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return 0;
    if (!point_ok(pt) || Q->input_points.size() >= wrapper::kMaxPoints)
        return Q->input_points.size();
    Q->insert_parallel(pt);
    return Q->input_points.size();
}

size_t qdt_insert_constraint(qdt_t* ptr, size_t v0, size_t v1)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return 0;
    return Q->insert_constraint(v0, v1) ? 1 : 0;
}

size_t qdt_insert_constraints_batch(qdt_t* ptr, const qdt_edge* edges, size_t num_edges)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return 0;
    return Q->insert_constraints_batch(edges, num_edges);
}

size_t qdt_triangles_size(const qdt_t* ptr)
{
    auto Q = reinterpret_cast<const wrapper::Triangulation*>(ptr);
    if (Q)
        return Q->triangles.size();
    return 0;
}

size_t qdt_get_triangle(const qdt_t* ptr, size_t T[3], size_t id)
{
    auto Q = reinterpret_cast<const wrapper::Triangulation*>(ptr);
    if (!Q || !T)
        return 0;
    if (id >= Q->triangles.size())
        return 0;
    for (int i = 0; i < 3; i++)
        T[i] = Q->triangles[id][i];
    return 1;
}

size_t qdt_points_size(const qdt_t* ptr)
{
    auto Q = reinterpret_cast<const wrapper::Triangulation*>(ptr);
    if (Q)
        return Q->input_points.size();
    return 0;
}

size_t qdt_get_point(const qdt_t* ptr, qdt_vec2d* pt, size_t id)
{
    auto Q = reinterpret_cast<const wrapper::Triangulation*>(ptr);
    if (!Q || !pt)
        return 0;
    if (id >= Q->input_points.size())
        return 0;
    *pt = Q->input_points[id];
    return 1;
}

size_t qdt_point_is_used(const qdt_t* ptr, size_t id)
{
    auto Q = reinterpret_cast<const wrapper::Triangulation*>(ptr);
    if (!Q)
        return 0;
    if (id >= Q->input_points.size() || id >= Q->rep.size())
        return 0;
    // Without triangles nothing is referenced. With triangles, every distinct
    // point is a vertex of the triangulation, and of several points with the
    // same coordinates only the lowest index (the representative) is.
    if (Q->triangles.empty())
        return 0;
    return Q->rep[id] == id ? 1 : 0;
}

void qdt_destroy(qdt_t* ptr)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    delete Q;
}
