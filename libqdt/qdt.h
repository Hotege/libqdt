/*
 * Copyright (c) 2026 Hotege
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
 * Coordinates passed to QDT must be finite (no NaN or infinity).
 */
struct qdt_vec2d
{
    double x; /**< X coordinate. */
    double y; /**< Y coordinate. */
};
typedef struct qdt_vec2d qdt_vec2d;

/**
 * @brief Opaque handle to a Delaunay triangulation context.
 *
 * Create it with qdt_create(), run triangulation with qdt_execute() or
 * qdt_execute_parallel(), optionally add points with qdt_insert_point() or
 * qdt_insert_point_parallel(), query results with qdt_triangles_size() and
 * qdt_get_triangle(), and release it with qdt_destroy().
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
 * added with qdt_insert_point(). The points are copied, so pts is not
 * modified and does not need to stay valid after the call.
 *
 * Triangle vertex indices in the result refer to positions in pts. If several
 * points have identical coordinates, only the one with the lowest index
 * appears in the triangles. With fewer than three distinct points, or with
 * all points collinear, the result has no triangles.
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pts Array of n 2D points with finite coordinates.
 * @param n Number of points in pts.
 */
QDT_API void qdt_execute(qdt_t* ptr, qdt_vec2d* pts, size_t n);

/**
 * @brief Computes the Delaunay triangulation using a parallel implementation.
 *
 * Produces the same set of triangles as qdt_execute(), but not necessarily in
 * the same order. The number of threads is std::thread::hardware_concurrency()
 * capped at 64. With fewer than 200000 distinct points the computation runs
 * on the calling thread, and for moderately sized inputs it may be no faster
 * than qdt_execute().
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pts Array of n 2D points with finite coordinates.
 * @param n Number of points in pts.
 */
QDT_API void qdt_execute_parallel(qdt_t* ptr, qdt_vec2d* pts, size_t n);

/**
 * @brief Inserts a single 2D point into the triangulation context.
 *
 * The point is appended to the points of the context and the triangulation is
 * updated incrementally. Its index is the return value minus 1.
 *
 * If pt has the same coordinates as an existing point, it still receives an
 * index, but no triangle references it.
 *
 * The cost of one call is proportional to the current number of points, and
 * the first call after qdt_execute() additionally builds an internal mesh
 * (extra time and memory). Use this function to add a small number of points;
 * to change many points, call qdt_execute() again with the full point set.
 * Triangle ids obtained before the call are no longer valid afterwards.
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
 * is single-threaded.
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pt 2D point to insert, with finite coordinates.
 * @return The total number of points in the context after the call
 *         (including duplicates), or 0 if ptr is NULL.
 */
QDT_API size_t qdt_insert_point_parallel(qdt_t* ptr, qdt_vec2d pt);

/**
 * @brief Returns the current number of triangles in the context.
 *
 * Reflects the last qdt_execute(), qdt_execute_parallel(),
 * qdt_insert_point() or qdt_insert_point_parallel() call.
 *
 * @param ptr Triangulation context.
 * @return Number of triangles, or 0 if ptr is NULL. Valid triangle ids are in
 *         [0, qdt_triangles_size(ptr)).
 */
QDT_API size_t qdt_triangles_size(qdt_t* ptr);

/**
 * @brief Retrieves the vertex indices of one triangle.
 *
 * The vertices are listed counter-clockwise. Which vertex comes first is
 * unspecified. Triangle ids are not stable: they may change after any call to
 * qdt_execute(), qdt_execute_parallel(), qdt_insert_point() or
 * qdt_insert_point_parallel().
 *
 * @param ptr Triangulation context.
 * @param T Output array receiving three point indices. For the batch
 *          functions the indices refer to the point array passed to
 *          qdt_execute() or qdt_execute_parallel(). Points added by
 *          qdt_insert_point() or qdt_insert_point_parallel() are indexed in
 *          insertion order, following the existing points.
 * @param id Triangle id, in [0, qdt_triangles_size(ptr)).
 * @return 1 if the triangle exists and T was filled, 0 otherwise (T is not
 *         modified).
 */
QDT_API size_t qdt_get_triangle(qdt_t* ptr, size_t T[3], size_t id);

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