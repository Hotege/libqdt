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
 * qdt_execute_parallel(), query results with qdt_triangles_size() and
 * qdt_get_triangle(), and release it with qdt_destroy().
 */
typedef struct qdt_t qdt_t;

/**
 * @brief Creates a new triangulation context.
 *
 * @return A non-NULL handle on success, or NULL if memory allocation fails.
 */
QDT_API qdt_t* qdt_create();

/**
 * @brief Computes the Delaunay triangulation of a 2D point set.
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pts Array of n 2D points.
 * @param n Number of points in pts.
 */
QDT_API void qdt_execute(qdt_t* ptr, qdt_vec2d* pts, size_t n);

/**
 * @brief Computes the Delaunay triangulation using a parallel implementation.
 *
 * The result is the same as qdt_execute(), but the computation may use
 * multiple threads for better performance on large point sets.
 *
 * @param ptr Triangulation context created by qdt_create().
 * @param pts Array of n 2D points.
 * @param n Number of points in pts.
 */
QDT_API void qdt_execute_parallel(qdt_t* ptr, qdt_vec2d* pts, size_t n);

/**
 * @brief Returns the number of triangles produced by the last execution.
 *
 * @param ptr Triangulation context.
 * @return Number of triangles. Valid triangle ids are in
 *         [0, qdt_triangles_size(ptr)).
 */
QDT_API size_t qdt_triangles_size(qdt_t* ptr);

/**
 * @brief Retrieves the vertex indices of one triangle.
 *
 * @param ptr Triangulation context.
 * @param T Output array receiving three point indices. T[0], T[1], and T[2]
 *          are indices into the point array passed to qdt_execute() or
 *          qdt_execute_parallel().
 * @param id Triangle index, in [0, qdt_triangles_size(ptr)).
 */
QDT_API void qdt_get_triangle(qdt_t* ptr, size_t T[3], size_t id);

/**
 * @brief Destroys a triangulation context and releases all associated resources.
 *
 * @param ptr Context to destroy.
 */
QDT_API void qdt_destroy(qdt_t* ptr);

#if __cplusplus
}
#endif

#endif