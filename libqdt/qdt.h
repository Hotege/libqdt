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

/**
 * @brief The two end points of a constraint, as point indices.
 */
struct qdt_edge
{
    size_t v0; /**< Index of the first end point. */
    size_t v1; /**< Index of the second end point. */
};
typedef struct qdt_edge qdt_edge;

/**
 * @brief Result of a QDT function that can fail.
 *
 * QDT_OK is 0 and every error is negative, so `status != QDT_OK` (or
 * `status < 0`) tests for failure. qdt_status_string() gives a description.
 *
 * What an error does to the context:
 *  - The call changes nothing, and the triangulation is exactly as it was:
 *    QDT_ERR_NULL_CONTEXT from any function; QDT_ERR_NOT_FINITE and
 *    QDT_ERR_TOO_MANY_POINTS from qdt_insert_point() and
 *    qdt_insert_point_parallel(); QDT_ERR_NULL_ARGUMENT from
 *    qdt_insert_constraints_batch(); and every rejection of a constraint
 *    (QDT_ERR_INVALID_INDEX up to QDT_ERR_THROUGH_POINT).
 *  - The context is emptied (no points, no triangles, no constraints):
 *    QDT_ERR_TOO_MANY_POINTS, QDT_ERR_NULL_ARGUMENT and QDT_ERR_NOT_FINITE
 *    from qdt_execute() and qdt_execute_parallel(), which always replace the
 *    previous contents, and QDT_ERR_OUT_OF_MEMORY and QDT_ERR_INTERNAL from any
 *    function. A function that fails halfway may have changed some internal
 *    data and not other data, so the context is reset instead of being kept in
 *    an inconsistent state.
 *
 * An emptied context stays valid. It can be queried, reused (for example by
 * calling qdt_execute() again) and destroyed.
 */
typedef enum qdt_status
{
    QDT_OK = 0,                       /**< Success. */
    QDT_ERR_NULL_CONTEXT = -1,        /**< The context pointer is NULL. */
    QDT_ERR_NULL_ARGUMENT = -2,       /**< A required pointer argument is NULL: pts with n > 0, or edges with num_edges > 0. */
    QDT_ERR_TOO_MANY_POINTS = -3,     /**< More than 100000000 points. */
    QDT_ERR_NOT_FINITE = -4,          /**< A coordinate is NaN or infinite. */
    QDT_ERR_INVALID_INDEX = -5,       /**< A constraint end point is not a point of the context. */
    QDT_ERR_SAME_POINT = -6,          /**< The end points of a constraint are the same point (the same index, or points with identical coordinates). */
    QDT_ERR_NO_TRIANGLES = -7,        /**< A constraint was requested, but the context has no triangles (fewer than three distinct points, or all points collinear). */
    QDT_ERR_CROSSES_CONSTRAINT = -8,  /**< The constraint would cross an edge that is already constrained. */
    QDT_ERR_THROUGH_POINT = -9,       /**< The constraint segment passes exactly through a third point. */
    QDT_ERR_OUT_OF_MEMORY = -10,      /**< Memory ran out. The context is emptied. */
    QDT_ERR_INTERNAL = -11            /**< Unexpected internal failure. The context is emptied. */
} qdt_status;

/**
 * @brief Returns a short English description of a status code.
 *
 * @param status A value returned by one of the functions of this interface.
 * @return A string with static storage duration (never NULL). An unknown
 *         value gives "unknown status".
 */
QDT_API const char* qdt_status_string(qdt_status status);

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
 * Passing a NULL context to any function is safe. The functions that return a
 * qdt_status return QDT_ERR_NULL_CONTEXT; the query functions return 0.
 *
 * No function of this interface lets a C++ exception escape, so all of them
 * can be called from C. Failures, including running out of memory, are
 * reported through qdt_status; see there for what happens to the context.
 */
typedef struct qdt_t qdt_t;

/**
 * @brief Creates a new, empty triangulation context.
 *
 * @return A handle that must be released with qdt_destroy(), or NULL if the
 *         memory for the context cannot be allocated.
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
 * Invalid input is rejected with an error code. The checks run in this order:
 * QDT_ERR_TOO_MANY_POINTS (more than 100000000 points), QDT_ERR_NULL_ARGUMENT
 * (pts is NULL with n > 0) and QDT_ERR_NOT_FINITE (a NaN or infinite
 * coordinate in pts). A rejected call still replaces the previous contents, so
 * the context is empty afterwards (no points, no triangles). Running out of
 * memory gives QDT_ERR_OUT_OF_MEMORY and the same empty context.
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pts Array of n 2D points with finite coordinates.
 * @param n Number of points in pts.
 * @return QDT_OK, QDT_ERR_NULL_CONTEXT (nothing happens), or one of the error
 *         codes above.
 */
QDT_API qdt_status qdt_execute(qdt_t* ptr, const qdt_vec2d* pts, size_t n);

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
 * the serial algorithm and still returns the complete result with QDT_OK.
 * Invalid input and running out of memory are handled as in qdt_execute().
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pts Array of n 2D points with finite coordinates.
 * @param n Number of points in pts.
 * @return The same values as qdt_execute().
 */
QDT_API qdt_status qdt_execute_parallel(qdt_t* ptr, const qdt_vec2d* pts, size_t n);

/**
 * @brief Inserts a single 2D point into the triangulation context.
 *
 * The point is appended to the points of the context and the triangulation is
 * updated incrementally. Its index is qdt_points_size(ptr) - 1 after the call;
 * pass an index pointer to receive it directly.
 *
 * If pt has the same coordinates as an existing point, it still receives an
 * index, but no triangle references it.
 *
 * A point with a NaN or infinite coordinate (QDT_ERR_NOT_FINITE), or one that
 * would take the context beyond 100000000 points (QDT_ERR_TOO_MANY_POINTS), is
 * rejected: nothing is changed.
 *
 * If memory runs out during the insertion, the call returns
 * QDT_ERR_OUT_OF_MEMORY and the context is emptied (see qdt_status).
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
 * @param index Optional, may be NULL. On success it receives the index of the
 *              new point (duplicates included). It is not written when the
 *              call fails.
 * @return QDT_OK, QDT_ERR_NULL_CONTEXT, QDT_ERR_NOT_FINITE,
 *         QDT_ERR_TOO_MANY_POINTS, QDT_ERR_OUT_OF_MEMORY or QDT_ERR_INTERNAL.
 */
QDT_API qdt_status qdt_insert_point(qdt_t* ptr, qdt_vec2d pt, size_t* index);

/**
 * @brief Inserts a single 2D point, building the internal mesh in parallel.
 *
 * Same effect and result as qdt_insert_point(). Only the one-time mesh build
 * performed by the first insertion after a batch call uses multiple threads
 * (when the triangulation has at least 50000 triangles). The insertion itself
 * is single-threaded. A rejected point and running out of memory are handled
 * as in qdt_insert_point().
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pt 2D point to insert, with finite coordinates.
 * @param index Optional, may be NULL; see qdt_insert_point().
 * @return The same values as qdt_insert_point().
 */
QDT_API qdt_status qdt_insert_point_parallel(qdt_t* ptr, qdt_vec2d pt, size_t* index);

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
 * A constraint that cannot be inserted is rejected, the triangulation is left
 * exactly as it was, and the return value says why. The reasons are checked in
 * this order:
 *  - QDT_ERR_INVALID_INDEX: v0 or v1 is not a point of the context,
 *  - QDT_ERR_SAME_POINT: v0 and v1 are the same point,
 *  - QDT_ERR_NO_TRIANGLES: the context has no triangles (fewer than three
 *    distinct points, or all points collinear),
 *  - QDT_ERR_CROSSES_CONSTRAINT: the segment would cross an edge that is
 *    already constrained,
 *  - QDT_ERR_THROUGH_POINT: the segment passes exactly through a third point
 *    (insert the pieces between consecutive points instead).
 * If the segment both crosses a constrained edge and passes through a point,
 * either of the last two may be reported.
 *
 * If memory runs out while the constraint is inserted, the call returns
 * QDT_ERR_OUT_OF_MEMORY, but unlike a rejection it does not leave the
 * triangulation as it was: the context is emptied (see qdt_status). The same
 * holds for QDT_ERR_INTERNAL.
 *
 * Constraints are kept until the next qdt_execute() or qdt_execute_parallel(),
 * which starts from scratch. Like qdt_insert_point(), the first call after a
 * batch call builds an internal mesh. Triangle ids obtained before the call
 * are no longer valid afterwards.
 *
 * @param ptr Triangulation context.
 * @param v0 Index of the first end point.
 * @param v1 Index of the second end point.
 * @return QDT_OK if the edge is now part of the triangulation as a constrained
 *         edge (including when it already was an edge); QDT_ERR_NULL_CONTEXT;
 *         one of the rejection reasons above; QDT_ERR_OUT_OF_MEMORY or
 *         QDT_ERR_INTERNAL.
 */
QDT_API qdt_status qdt_insert_constraint(qdt_t* ptr, size_t v0, size_t v1);

/**
 * @brief Inserts several constraints, one after another.
 *
 * Same as calling qdt_insert_constraint() for each entry in order: a rejected
 * entry does not stop the others, and a later entry that crosses an earlier
 * accepted one is rejected.
 *
 * If memory runs out (QDT_ERR_OUT_OF_MEMORY) or an internal error occurs
 * (QDT_ERR_INTERNAL), processing stops and the context is emptied (see
 * qdt_status); no entry counts as accepted then.
 *
 * @param ptr Triangulation context.
 * @param edges Array of num_edges end point pairs.
 * @param num_edges Number of entries in edges. Zero is allowed and does
 *                  nothing.
 * @param results Optional, may be NULL. If it is not NULL it must point to an
 *                array of num_edges entries, and entry i receives the status
 *                that entry i of edges got. When the whole call fails
 *                (QDT_ERR_NULL_CONTEXT, QDT_ERR_NULL_ARGUMENT,
 *                QDT_ERR_OUT_OF_MEMORY, QDT_ERR_INTERNAL), every entry
 *                receives that code.
 * @return QDT_OK if every entry was accepted (or already existed as an edge),
 *         and also when num_edges is 0. QDT_ERR_NULL_CONTEXT, or
 *         QDT_ERR_NULL_ARGUMENT when edges is NULL and num_edges > 0, with
 *         nothing changed. QDT_ERR_OUT_OF_MEMORY or QDT_ERR_INTERNAL when the
 *         call failed as a whole. Otherwise the status of the first rejected
 *         entry; the entries after it were still processed.
 */
QDT_API qdt_status qdt_insert_constraints_batch(qdt_t* ptr, const qdt_edge* edges, size_t num_edges, qdt_status* results);

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
 *         modified; this includes ptr or T being NULL). Unlike the functions
 *         that return a qdt_status, the query functions return 1 for success.
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