
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
        };

        bool                     mesh_built_ = false;
        std::vector<MeshVertex>  mesh_vertices_;
        std::vector<MeshHalfEdge> mesh_halfedges_;
        std::vector<int>         free_halfedges_;
        int                      last_face_he_ = 0;

        void ensure_mesh_();
        int  alloc_halfedge_();
        void free_halfedge_(int h);
        int  locate_(const qdt_vec2d& p, int& hull_he);
        void insert_by_rebuild_(const qdt_vec2d& p);
        void rebuild_output_();
    };
}

#endif
