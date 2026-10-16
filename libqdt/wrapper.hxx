#if !defined(_LIB_QDT_WRAPPER_HXX_)
#define _LIB_QDT_WRAPPER_HXX_

#include <array>
#include <vector>

namespace wrapper
{
    struct Triangulation
    {
        std::vector<std::array<size_t, 3>> triangles;
        std::vector<size_t> rep;
        int unique_points = 0;

        std::vector<qdt_vec2d> input_points;

        void insert(const qdt_vec2d& p);
        void insert_parallel(const qdt_vec2d& p);

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
            int origin;
            int twin;
            int next;
            int triangle_idx;
        };

        bool                      mesh_built_ = false;
        std::vector<MeshVertex>   mesh_vertices_;
        std::vector<MeshHalfEdge> mesh_halfedges_;
        std::vector<int>          free_halfedges_;
        int                       last_face_he_ = 0;
        std::vector<int>          triangle_to_he_;

        void insert_after_mesh_(const qdt_vec2d& p);
        int  alloc_halfedge_();
        void free_halfedge_(int h);
        int  locate_(const qdt_vec2d& p, int& hull_he);
        void insert_by_rebuild_(const qdt_vec2d& p);
        size_t add_triangle_to_output_(size_t id0, size_t id1, size_t id2, int h0);
        void remove_triangle_by_index_(size_t idx);
    };
}

#endif