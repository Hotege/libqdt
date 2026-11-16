#include <cmath>
#include <new>

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

    // Checks the input of a batch call: at most wrapper::kMaxPoints points, a
    // non-NULL array unless it is empty, and finite coordinates only. Returns
    // QDT_OK or the first violation found, in that order.
    qdt_status batch_input_status(const qdt_vec2d* pts, size_t n)
    {
        if (n > wrapper::kMaxPoints)
            return QDT_ERR_TOO_MANY_POINTS;
        if (n == 0)
            return QDT_OK;
        if (!pts)
            return QDT_ERR_NULL_ARGUMENT;
        for (size_t i = 0; i < n; ++i)
        {
            if (!point_ok(pts[i]))
                return QDT_ERR_NOT_FINITE;
        }
        return QDT_OK;
    }

    // Puts the context into the empty state: no points, no triangles, no
    // constraints and no internal mesh.
    //
    // Every function below that changes a context runs inside try/catch so
    // that no C++ exception reaches the C caller. A failed operation (in
    // practice std::bad_alloc) may already have updated some of the context's
    // containers but not others, so what is left is not a usable
    // triangulation. The context is therefore reset instead of being kept in
    // a half-changed state. Move-assigning a default-constructed object does
    // not allocate and cannot throw, which makes this safe to call from a
    // catch handler.
    void reset_context(wrapper::Triangulation* Q) noexcept
    {
        *Q = wrapper::Triangulation();
    }

    // Maps the exception that is being handled to a status. Call it only from
    // inside a catch block.
    qdt_status status_from_current_exception() noexcept
    {
        try
        {
            throw;
        }
        catch (const std::bad_alloc&)
        {
            return QDT_ERR_OUT_OF_MEMORY;
        }
        catch (...)
        {
            return QDT_ERR_INTERNAL;
        }
    }

    void fill_results(qdt_status* results, size_t n, qdt_status status) noexcept
    {
        if (!results)
            return;
        for (size_t i = 0; i < n; ++i)
            results[i] = status;
    }

    // Shared by qdt_execute() and qdt_execute_parallel().
    qdt_status execute_impl(qdt_t* ptr, const qdt_vec2d* pts, size_t n, bool parallel_run)
    {
        auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
        if (!Q)
            return QDT_ERR_NULL_CONTEXT;
        try
        {
            const qdt_status input = batch_input_status(pts, n);
            if (input != QDT_OK)
            {
                // Rejected input: the call still replaces the previous contents,
                // so the context is left empty rather than holding stale results.
                reset_context(Q);
                return input;
            }
            *Q = parallel_run ? parallel::delaunay_triangulate(pts, n)
                : wrapper::delaunay_triangulate(pts, n);
            return QDT_OK;
        }
        catch (...)
        {
            const qdt_status st = status_from_current_exception();
            reset_context(Q);
            return st;
        }
    }

    // Shared by qdt_insert_point() and qdt_insert_point_parallel().
    qdt_status insert_point_impl(qdt_t* ptr, const qdt_vec2d& pt, size_t* index, bool parallel_mesh)
    {
        auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
        if (!Q)
            return QDT_ERR_NULL_CONTEXT;
        try
        {
            // A rejected point is not added and nothing is changed.
            if (!point_ok(pt))
                return QDT_ERR_NOT_FINITE;
            if (Q->input_points.size() >= wrapper::kMaxPoints)
                return QDT_ERR_TOO_MANY_POINTS;
            if (parallel_mesh)
                Q->insert_parallel(pt);
            else
                Q->insert(pt);
            if (index)
                *index = Q->input_points.size() - 1;
            return QDT_OK;
        }
        catch (...)
        {
            const qdt_status st = status_from_current_exception();
            reset_context(Q);
            return st;
        }
    }
}

const char* qdt_status_string(qdt_status status)
{
    switch (status)
    {
    case QDT_OK:                     return "success";
    case QDT_ERR_NULL_CONTEXT:       return "the context is NULL";
    case QDT_ERR_NULL_ARGUMENT:      return "a required pointer argument is NULL";
    case QDT_ERR_TOO_MANY_POINTS:    return "more than 100000000 points";
    case QDT_ERR_NOT_FINITE:         return "a coordinate is NaN or infinite";
    case QDT_ERR_INVALID_INDEX:      return "an end point index is not a point of the context";
    case QDT_ERR_SAME_POINT:         return "both end points are the same point";
    case QDT_ERR_NO_TRIANGLES:       return "the context has no triangles";
    case QDT_ERR_CROSSES_CONSTRAINT: return "the segment crosses an existing constraint";
    case QDT_ERR_THROUGH_POINT:      return "the segment passes through another point";
    case QDT_ERR_OUT_OF_MEMORY:      return "out of memory";
    case QDT_ERR_INTERNAL:           return "internal error";
    }
    return "unknown status";
}

qdt_t* qdt_create()
{
    // A default-constructed Triangulation does not allocate, so the only
    // allocation that can fail here is the object itself. nothrow new turns
    // that failure into NULL instead of std::bad_alloc.
    return reinterpret_cast<qdt_t*>(new (std::nothrow) wrapper::Triangulation);
}

qdt_status qdt_execute(qdt_t* ptr, const qdt_vec2d* pts, size_t n)
{
    return execute_impl(ptr, pts, n, false);
}

qdt_status qdt_execute_parallel(qdt_t* ptr, const qdt_vec2d* pts, size_t n)
{
    return execute_impl(ptr, pts, n, true);
}

qdt_status qdt_insert_point(qdt_t* ptr, qdt_vec2d pt, size_t* index)
{
    return insert_point_impl(ptr, pt, index, false);
}

qdt_status qdt_insert_point_parallel(qdt_t* ptr, qdt_vec2d pt, size_t* index)
{
    return insert_point_impl(ptr, pt, index, true);
}

qdt_status qdt_insert_constraint(qdt_t* ptr, size_t v0, size_t v1)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
        return QDT_ERR_NULL_CONTEXT;
    try
    {
        const qdt_status st = Q->insert_constraint(v0, v1);
        // An internal inconsistency may have left the triangulation changed.
        if (st == QDT_ERR_INTERNAL)
            reset_context(Q);
        return st;
    }
    catch (...)
    {
        const qdt_status st = status_from_current_exception();
        reset_context(Q);
        return st;
    }
}

qdt_status qdt_insert_constraints_batch(qdt_t* ptr, const qdt_edge* edges, size_t num_edges, qdt_status* results)
{
    auto Q = reinterpret_cast<wrapper::Triangulation*>(ptr);
    if (!Q)
    {
        fill_results(results, num_edges, QDT_ERR_NULL_CONTEXT);
        return QDT_ERR_NULL_CONTEXT;
    }
    if (num_edges == 0)
        return QDT_OK;
    if (!edges)
    {
        fill_results(results, num_edges, QDT_ERR_NULL_ARGUMENT);
        return QDT_ERR_NULL_ARGUMENT;
    }
    try
    {
        const qdt_status st = Q->insert_constraints_batch(edges, num_edges, results);
        if (st == QDT_ERR_INTERNAL)
            reset_context(Q);
        return st;
    }
    catch (...)
    {
        // The context is empty now, so no entry counts as accepted.
        const qdt_status st = status_from_current_exception();
        reset_context(Q);
        fill_results(results, num_edges, st);
        return st;
    }
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