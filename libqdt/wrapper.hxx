#if !defined(_LIB_QDT_WRAPPER_HXX_)
#define _LIB_QDT_WRAPPER_HXX_

#include <array>
#include <vector>
#include <unordered_set>
#include "qdt.h"

namespace wrapper
{
    // Largest number of points a context accepts. The internal indices are
    // 32-bit: the parallel builder allocates 5 * n + 64 Ki quad-edges and a
    // quad-edge id is 4 * k + r, which stays below 2^31 up to about 107 million
    // points. Larger inputs are rejected at the API boundary (qdt.cpp).
    constexpr size_t kMaxPoints = 100000000;

    struct Triangulation
    {
        std::vector<std::array<size_t, 3>> triangles;
        std::vector<size_t> rep;
        int unique_points = 0;

        std::vector<qdt_vec2d> input_points;

        void insert(const qdt_vec2d& p);
        void insert_parallel(const qdt_vec2d& p);

        bool insert_constraint(size_t orig_v0, size_t orig_v1);

        size_t insert_constraints_batch(const qdt_edge* edges, size_t num_edges);

        void ensure_mesh_();
        void ensure_mesh_parallel_();

    private:
        struct MeshVertex
        {
            qdt_vec2d p;
            size_t    orig_id;
            int       halfedge;
        };

        struct MeshHalfEdge
        {
            int  origin;
            int  twin;
            int  next;
            int  triangle_idx;
            bool is_constrained = false;
        };

        bool                      mesh_built_ = false;
        std::vector<MeshVertex>   mesh_vertices_;
        // Index of the mesh vertex for each original point id (-1: the id has
        // no mesh vertex, e.g. a duplicate of another point). Built together
        // with the mesh and extended by every insertion; it may be shorter than
        // input_points when the latest points are duplicates. Mesh vertices are
        // only ever appended, so an entry stays valid until the mesh is rebuilt.
        std::vector<int>          orig_to_local_;
        std::vector<MeshHalfEdge> mesh_halfedges_;
        std::vector<int>          free_halfedges_;
        int                       last_face_he_ = 0;
        std::vector<int>          triangle_to_he_;

        void insert_after_mesh_(const qdt_vec2d& p);
        int  alloc_halfedge_();
        void free_halfedge_(int h);
        int  locate_(const qdt_vec2d& p, int& hull_he);
        int  locate_scan_(const qdt_vec2d& p, int& hull_he);
        void insert_by_rebuild_(const qdt_vec2d& p);
        size_t add_triangle_to_output_(size_t id0, size_t id1, size_t id2, int h0);
        void remove_triangle_by_index_(size_t idx);

        template <class F>
        void for_each_out_(int x_local, F&& f) const;
        int  find_halfedge_(int u_local, int v_local) const;
        bool flip_edge_(int h);
        bool segment_intersects_(const qdt_vec2d& a, const qdt_vec2d& b, const qdt_vec2d& c, const qdt_vec2d& d) const;
        int  orig_to_local_vtx_(size_t orig_id) const;
        void mark_constrained_(int h);
        bool trace_segment_(int u_local, int v_local, std::vector<int>& crossed) const;
        void legalize_(std::vector<int>& stack);
        bool insert_constraint_local_(int u_local, int v_local);
    };

    // Serial triangulation; defined in GuibasStolfiDwyer_Wrapper.cpp. Declared
    // here so that the parallel implementation can fall back to it.
    Triangulation delaunay_triangulate(const qdt_vec2d* pts, size_t sz);
}

#endif