/*
 * Copyright (c) 2026 QDT contributors
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the QDT library.
 */
#if !defined(_LIB_QDT_H_)
#define _LIB_QDT_H_

#include <stddef.h>

#if defined(_WIN32)
#if defined(QDT_DLL)
#define QDT_API __declspec(dllexport)
#else
#define QDT_API __declspec(dllimport)
#endif
#else
#define QDT_API
#endif

#if __cplusplus
extern "C"
{
#endif

/**
 * @brief A 2D point/vector with double-precision coordinates.
 *
 * Coordinates passed to QDT must be finite (no NaN or infinity). The
 * geometric tests are exact for any overall scale of the input (a point set
 * scaled to 1e-300 or to 1e300 works like one at unit size) and for sets that
 * mix very different magnitudes. The exception is subnormal numbers: nonzero
 * coordinates with an absolute value below about 2.2e-308 have too few
 * significant bits, and the output may then fail to be Delaunay without any
 * error being reported.
 */
struct qdt_vec2d
{
    double x; /**< X coordinate. */
    double y; /**< Y coordinate. */
};
typedef struct qdt_vec2d qdt_vec2d;

struct qdt_edge
{
    size_t v0;
    size_t v1;
};
typedef struct qdt_edge qdt_edge;

/**
 * @brief Opaque handle to a Delaunay triangulation context.
 *
 * Create it with qdt_create(), run triangulation with qdt_execute() or
 * qdt_execute_parallel(), optionally add points with qdt_insert_point() or
 * qdt_insert_point_parallel(), force edges into the triangulation with
 * qdt_insert_constraint() or qdt_insert_constraints_batch(), query results with
 * qdt_triangles_size() and qdt_get_triangle() (and the points with
 * qdt_points_size(), qdt_get_point() and qdt_point_is_used()), and release it
 * with qdt_destroy().
 *
 * A single context must not be used from several threads at the same time.
 * Different contexts are independent and may be used concurrently.
 *
 * Passing a NULL context to any function is safe: the function does nothing
 * and returns 0.
 */
typedef struct qdt_t qdt_t;

/**
 * @brief Creates a new, empty triangulation context.
 *
 * @note QDT does not catch C++ exceptions at this interface. Allocation
 *       failure is reported by throwing std::bad_alloc, not by returning
 *       NULL. The same applies to qdt_execute(), qdt_execute_parallel(),
 *       qdt_insert_point() and qdt_insert_point_parallel(). Letting an
 *       exception propagate into C code is undefined behavior.
 *
 * @return A handle that must be released with qdt_destroy().
 */
QDT_API qdt_t* qdt_create(void);

/**
 * @brief Computes the Delaunay triangulation of a 2D point set.
 *
 * Replaces everything previously stored in the context, including points
 * added with qdt_insert_point() and constraints. The points are copied, so pts is not
 * modified and does not need to stay valid after the call.
 *
 * Triangle vertex indices in the result refer to positions in pts. If several
 * points have identical coordinates, only the one with the lowest index
 * appears in the triangles. With fewer than three distinct points, or with
 * all points collinear, the result has no triangles.
 *
 * Invalid input is rejected: more than 100000000 points, a NULL pts with
 * n > 0, or a NaN or infinite coordinate. A rejected call still replaces the
 * previous contents, so the context is empty afterwards (no points, no
 * triangles) and an error is not reported in any other way.
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pts Array of n 2D points with finite coordinates.
 * @param n Number of points in pts.
 */
QDT_API void qdt_execute(qdt_t* ptr, const qdt_vec2d* pts, size_t n);

/**
 * @brief Computes the Delaunay triangulation using a parallel implementation.
 *
 * Produces the same set of triangles as qdt_execute(), but not necessarily in
 * the same order. The number of threads is std::thread::hardware_concurrency()
 * capped at 64. With fewer than 200000 distinct points the computation runs
 * on the calling thread, and for moderately sized inputs it may be no faster
 * than qdt_execute().
 *
 * If the parallel implementation cannot run (its pre-allocated edge storage is
 * exhausted, or the worker threads cannot be started), the call falls back to
 * the serial algorithm and still returns the complete result. Invalid input is
 * handled as in qdt_execute().
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pts Array of n 2D points with finite coordinates.
 * @param n Number of points in pts.
 */
QDT_API void qdt_execute_parallel(qdt_t* ptr, const qdt_vec2d* pts, size_t n);

/**
 * @brief Inserts a single 2D point into the triangulation context.
 *
 * The point is appended to the points of the context and the triangulation is
 * updated incrementally. Its index is the return value minus 1.
 *
 * If pt has the same coordinates as an existing point, it still receives an
 * index, but no triangle references it.
 *
 * A point with a NaN or infinite coordinate, or one that would take the
 * context beyond 100000000 points, is rejected: nothing is changed and the
 * current number of points is returned, so a return value that did not grow
 * signals the rejection.
 *
 * The cost of one call is proportional to the current number of points, and
 * the first call after qdt_execute() additionally builds an internal mesh
 * (extra time and memory). Use this function to add a small number of points;
 * to change many points, call qdt_execute() again with the full point set.
 * Triangle ids obtained before the call are no longer valid afterwards.
 *
 * Constraints (see qdt_insert_constraint()) are kept. If pt lies exactly on a
 * constrained edge, that edge is split at pt and both halves stay constrained.
 * The result is the constrained Delaunay triangulation of all points and
 * constraints.
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pt 2D point to insert, with finite coordinates.
 * @return The total number of points in the context after the call
 *         (including duplicates), or 0 if ptr is NULL.
 */
QDT_API size_t qdt_insert_point(qdt_t* ptr, qdt_vec2d pt);

/**
 * @brief Inserts a single 2D point, building the internal mesh in parallel.
 *
 * Same effect and result as qdt_insert_point(). Only the one-time mesh build
 * performed by the first insertion after a batch call uses multiple threads
 * (when the triangulation has at least 50000 triangles). The insertion itself
 * is single-threaded. A rejected point is handled as in qdt_insert_point().
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pt 2D point to insert, with finite coordinates.
 * @return The total number of points in the context after the call
 *         (including duplicates), or 0 if ptr is NULL.
 */
QDT_API size_t qdt_insert_point_parallel(qdt_t* ptr, qdt_vec2d pt);

/**
 * @brief Forces the segment between two points to be an edge of the
 *        triangulation.
 *
 * v0 and v1 are point indices, as in the triangles (an index of a duplicate
 * point acts as the lowest index with the same coordinates). The edge is made
 * a constrained edge: afterwards it is part of the output, no other edge
 * crosses it, and qdt_insert_point() never removes it. All other edges stay
 * locally Delaunay, so the result is the constrained Delaunay triangulation of
 * the points and all constraints inserted so far. With points in general
 * position that triangulation is unique, so the order in which constraints are
 * inserted does not matter; with cocircular points, equally valid results can
 * differ. The order of v0 and v1 does not matter either (also for edges on the
 * convex hull).
 *
 * The constraint is rejected, and the triangulation is left exactly as it was,
 * when
 *  - it would cross an edge that is already constrained,
 *  - the segment passes exactly through a third point (insert the pieces
 *    between consecutive points instead),
 *  - v0 and v1 are the same point, or an index is not a point of the context,
 *  - the context has no triangles (fewer than three distinct points, or all
 *    points collinear).
 *
 * Constraints are kept until the next qdt_execute() or qdt_execute_parallel(),
 * which starts from scratch. Like qdt_insert_point(), the first call after a
 * batch call builds an internal mesh. Triangle ids obtained before the call
 * are no longer valid afterwards.
 *
 * @param ptr Triangulation context.
 * @param v0 Index of the first end point.
 * @param v1 Index of the second end point.
 * @return 1 if the edge is now part of the triangulation as a constrained
 *         edge (including when it already was an edge), 0 if it was rejected
 *         or ptr is NULL.
 */
QDT_API size_t qdt_insert_constraint(qdt_t* ptr, size_t v0, size_t v1);

/**
 * @brief Inserts several constraints, one after another.
 *
 * Same as calling qdt_insert_constraint() for each entry in order: a rejected
 * entry does not stop the others, and a later entry that crosses an earlier
 * accepted one is rejected.
 *
 * @param ptr Triangulation context.
 * @param edges Array of num_edges end point pairs.
 * @param num_edges Number of entries in edges.
 * @return The number of entries that were accepted (or already existed as
 *         edges). 0 if ptr or edges is NULL.
 */
QDT_API size_t qdt_insert_constraints_batch(qdt_t* ptr, const qdt_edge* edges, size_t num_edges);

/**
 * @brief Returns the current number of triangles in the context.
 *
 * Reflects the last qdt_execute(), qdt_execute_parallel(),
 * qdt_insert_point(), qdt_insert_point_parallel() or constraint call.
 *
 * @param ptr Triangulation context.
 * @return Number of triangles, or 0 if ptr is NULL. Valid triangle ids are in
 *         [0, qdt_triangles_size(ptr)).
 */
QDT_API size_t qdt_triangles_size(const qdt_t* ptr);

/**
 * @brief Retrieves the vertex indices of one triangle.
 *
 * The vertices are listed counter-clockwise. Which vertex comes first is
 * unspecified. Triangle ids are not stable: they may change after any call to
 * qdt_execute(), qdt_execute_parallel(), qdt_insert_point(),
 * qdt_insert_point_parallel(), qdt_insert_constraint() or
 * qdt_insert_constraints_batch().
 *
 * @param ptr Triangulation context.
 * @param T Output array receiving three point indices. For the batch
 *          functions the indices refer to the point array passed to
 *          qdt_execute() or qdt_execute_parallel(). Points added by
 *          qdt_insert_point() or qdt_insert_point_parallel() are indexed in
 *          insertion order, following the existing points.
 * @param id Triangle id, in [0, qdt_triangles_size(ptr)).
 * @return 1 if the triangle exists and T was filled, 0 otherwise (T is not
 *         modified; this includes ptr or T being NULL).
 */
QDT_API size_t qdt_get_triangle(const qdt_t* ptr, size_t T[3], size_t id);

/**
 * @brief Returns the number of points in the context.
 *
 * This is the number of points given to qdt_execute() or
 * qdt_execute_parallel() plus the number of points accepted by
 * qdt_insert_point() or qdt_insert_point_parallel() since then, duplicates
 * included. Point indices, as used in the triangles and in the constraint
 * functions, are in [0, qdt_points_size(ptr)).
 *
 * @param ptr Triangulation context.
 * @return Number of points, or 0 if ptr is NULL.
 */
QDT_API size_t qdt_points_size(const qdt_t* ptr);

/**
 * @brief Retrieves the coordinates of one point.
 *
 * Works for every index in [0, qdt_points_size(ptr)), also for duplicates and
 * for points that no triangle references.
 *
 * @param ptr Triangulation context.
 * @param pt Output receiving the coordinates.
 * @param id Point index.
 * @return 1 if the point exists and pt was filled, 0 otherwise (pt is not
 *         modified; this includes ptr or pt being NULL).
 */
QDT_API size_t qdt_get_point(const qdt_t* ptr, qdt_vec2d* pt, size_t id);

/**
 * @brief Tells whether a point is a vertex of the triangulation.
 *
 * A point is not referenced by any triangle when
 *  - it duplicates the coordinates of a point with a lower index (only the
 *    lowest index is used), or
 *  - the context has no triangles (fewer than three distinct points, or all
 *    points collinear).
 * Every other point is a vertex of at least one triangle. Constraints do not
 * change this.
 *
 * @param ptr Triangulation context.
 * @param id Point index.
 * @return 1 if at least one triangle references the point, 0 if not, if id is
 *         not a point of the context, or if ptr is NULL.
 */
QDT_API size_t qdt_point_is_used(const qdt_t* ptr, size_t id);

/**
 * @brief Destroys a triangulation context and releases all associated resources.
 *
 * @param ptr Context to destroy. NULL is allowed and does nothing.
 */
QDT_API void qdt_destroy(qdt_t* ptr);

#if __cplusplus
}
#endif

#endif