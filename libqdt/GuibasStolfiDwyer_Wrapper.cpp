#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <thread>
#include <unordered_map>
#include <vector>

#include "qdt.h"
#include "builder.hxx"
#include "robust.hxx"
#include "wrapper.hxx"

namespace wrapper
{
    namespace detail
    {
        // Edge storage of the serial builder: one vector that grows as needed,
        // and a free list of the records that deleteEdge() gave back. Nothing
        // runs concurrently. See builder.hxx for what a policy has to provide.
        class SerialPolicy
        {
        public:
            void init(std::vector<QE>& qe, size_t nverts) { qe.reserve(3 * nverts + 16); }

            int alloc(std::vector<QE>& qe)
            {
                if (free_head_ >= 0)
                {
                    const int k = free_head_;
                    free_head_ = qe[k].next[0];
                    return k;
                }
                const int k = static_cast<int>(qe.size());
                qe.emplace_back();
                return k;
            }

            void release(std::vector<QE>& qe, int k)
            {
                qe[k].vtx[0] = -1;
                qe[k].next[0] = free_head_;
                free_head_ = k;
            }

            int edge_count(const std::vector<QE>& qe) const { return static_cast<int>(qe.size()); }

            bool split_build(int, int) const { return false; }
            bool split_merge(int) const { return false; }

            template <class A, class B>
            auto fork_join(bool, A&& first, B&& second) -> std::pair<decltype(first()), decltype(second())>
            {
                auto x = first();
                auto y = second();
                return std::pair<decltype(first()), decltype(second())>(x, y);
            }

        private:
            int free_head_ = -1;
        };
    }

    size_t Triangulation::add_triangle_to_output_(size_t id0, size_t id1, size_t id2, int h0)
    {
        size_t new_idx = triangles.size();
        triangles.push_back({ id0, id1, id2 });
        triangle_to_he_.push_back(h0);
        return new_idx;
    }

    void Triangulation::remove_triangle_by_index_(size_t idx)
    {
        size_t last_idx = triangles.size() - 1;
        if (idx != last_idx)
        {
            triangles[idx] = triangles[last_idx];
            int moved_he = triangle_to_he_[last_idx];
            triangle_to_he_[idx] = moved_he;

            int cur = moved_he;
            do {
                mesh_halfedges_[cur].triangle_idx = static_cast<int>(idx);
                cur = mesh_halfedges_[cur].next;
            } while (cur != moved_he);
        }
        triangles.pop_back();
        triangle_to_he_.pop_back();
    }

    void Triangulation::ensure_mesh_()
    {
        if (mesh_built_) return;

        mesh_vertices_.clear();
        mesh_halfedges_.clear();
        free_halfedges_.clear();
        last_face_he_ = 0;
        triangle_to_he_.clear();

        const size_t nTri = triangles.size();
        const size_t nHE = 3 * nTri;

        std::vector<int> orig_to_local(input_points.size(), -1);
        mesh_vertices_.reserve(nTri * 3);

        for (size_t t = 0; t < nTri; ++t)
        {
            const auto& tri = triangles[t];
            for (int i = 0; i < 3; ++i)
            {
                const size_t orig = tri[i];
                int& slot = orig_to_local[orig];
                if (slot < 0)
                {
                    slot = static_cast<int>(mesh_vertices_.size());
                    MeshVertex v;
                    v.p = input_points[orig];
                    v.orig_id = orig;
                    v.halfedge = -1;
                    mesh_vertices_.push_back(v);
                }
            }
        }

        mesh_halfedges_.resize(nHE);
        triangle_to_he_.resize(nTri);

        MeshHalfEdge* __restrict he_ptr = mesh_halfedges_.data();
        for (size_t t = 0; t < nTri; ++t)
        {
            const auto& tri = triangles[t];
            const int v0 = orig_to_local[tri[0]];
            const int v1 = orig_to_local[tri[1]];
            const int v2 = orig_to_local[tri[2]];

            const int base = static_cast<int>(3 * t);

            he_ptr[base + 0] = { v0, -1, base + 1, static_cast<int>(t), false };
            he_ptr[base + 1] = { v1, -1, base + 2, static_cast<int>(t), false };
            he_ptr[base + 2] = { v2, -1, base + 0, static_cast<int>(t), false };

            if (mesh_vertices_[v0].halfedge < 0) mesh_vertices_[v0].halfedge = base + 0;
            if (mesh_vertices_[v1].halfedge < 0) mesh_vertices_[v1].halfedge = base + 1;
            if (mesh_vertices_[v2].halfedge < 0) mesh_vertices_[v2].halfedge = base + 2;

            triangle_to_he_[t] = base;
        }

        const int nV = static_cast<int>(mesh_vertices_.size());

        std::vector<int> offset(nV + 1, 0);
        for (size_t i = 0; i < nHE; ++i)
            ++offset[he_ptr[i].origin + 1];
        for (int v = 0; v < nV; ++v)
            offset[v + 1] += offset[v];

        std::vector<int> adj(nHE);
        {
            std::vector<int> cursor(offset.begin(), offset.begin() + nV);
            for (size_t i = 0; i < nHE; ++i)
            {
                const int u = he_ptr[i].origin;
                adj[cursor[u]++] = static_cast<int>(i);
            }
        }

        for (size_t i = 0; i < nHE; ++i)
        {
            if (he_ptr[i].twin >= 0) continue;

            const int u = he_ptr[i].origin;
            const int v = he_ptr[he_ptr[i].next].origin;

            const int start = offset[v];
            const int end = offset[v + 1];

            for (int k = start; k < end; ++k)
            {
                const int h2 = adj[k];
                if (he_ptr[h2].origin != v) continue;
                const int u2 = he_ptr[he_ptr[h2].next].origin;
                if (u2 == u)
                {
                    he_ptr[i].twin = h2;
                    he_ptr[h2].twin = static_cast<int>(i);
                    break;
                }
            }
        }

        orig_to_local_ = std::move(orig_to_local);
        mesh_built_ = true;
    }

    namespace
    {
        // Joins every thread of a vector when it goes out of scope.
        //
        // Starting a thread allocates, so emplace_back() on the workers can
        // throw after some of them are already running. Destroying a joinable
        // std::thread calls std::terminate(), so without this guard such an
        // exception would end the process instead of reaching the catch
        // handler at the C API boundary. Declare it right after the vector, and
        // declare everything the threads refer to by reference before it: it is
        // then destroyed first during stack unwinding, so the threads have
        // finished before those objects are destroyed.
        class ThreadJoiner
        {
        public:
            explicit ThreadJoiner(std::vector<std::thread>& threads) : threads_(threads) {}
            ~ThreadJoiner()
            {
                for (auto& t : threads_)
                    if (t.joinable()) t.join();
            }
            ThreadJoiner(const ThreadJoiner&) = delete;
            ThreadJoiner& operator=(const ThreadJoiner&) = delete;

        private:
            std::vector<std::thread>& threads_;
        };
    }

    void Triangulation::ensure_mesh_parallel_()
    {
        if (mesh_built_) return;

        unsigned nthreads = std::thread::hardware_concurrency();
        if (nthreads == 0) nthreads = 1;
        if (nthreads > 64) nthreads = 64;

        const size_t nTri = triangles.size();

        if (nthreads == 1 || nTri < 50000)
        {
            ensure_mesh_();
            return;
        }

        mesh_vertices_.clear();
        mesh_halfedges_.clear();
        free_halfedges_.clear();
        last_face_he_ = 0;
        triangle_to_he_.clear();

        const size_t nHE = 3 * nTri;

        std::vector<int> orig_to_local(input_points.size(), -1);
        mesh_vertices_.reserve(nTri * 3);

        for (size_t t = 0; t < nTri; ++t)
        {
            const auto& tri = triangles[t];
            for (int i = 0; i < 3; ++i)
            {
                const size_t orig = tri[i];
                int& slot = orig_to_local[orig];
                if (slot < 0)
                {
                    slot = static_cast<int>(mesh_vertices_.size());
                    MeshVertex v;
                    v.p = input_points[orig];
                    v.orig_id = orig;
                    v.halfedge = -1;
                    mesh_vertices_.push_back(v);
                }
            }
        }

        mesh_halfedges_.resize(nHE);
        triangle_to_he_.resize(nTri);

        const int nV = static_cast<int>(mesh_vertices_.size());

        std::unique_ptr<std::atomic<int>[]> v_he(new std::atomic<int>[nV]);
        for (int i = 0; i < nV; ++i)
            v_he[i].store(-1, std::memory_order_relaxed);

        {
            std::vector<std::thread> workers;
            workers.reserve(nthreads);
            ThreadJoiner joiner(workers);

            auto fill = [&](size_t t0, size_t t1) {
                MeshHalfEdge* hp = mesh_halfedges_.data();
                for (size_t t = t0; t < t1; ++t)
                {
                    const auto& tri = triangles[t];
                    const int v0 = orig_to_local[tri[0]];
                    const int v1 = orig_to_local[tri[1]];
                    const int v2 = orig_to_local[tri[2]];
                    const int base = static_cast<int>(3 * t);

                    hp[base + 0] = { v0, -1, base + 1, static_cast<int>(t), false };
                    hp[base + 1] = { v1, -1, base + 2, static_cast<int>(t), false };
                    hp[base + 2] = { v2, -1, base + 0, static_cast<int>(t), false };

                    triangle_to_he_[t] = base;

                    int expected = -1;
                    v_he[v0].compare_exchange_strong(
                        expected, base + 0, std::memory_order_relaxed);
                    expected = -1;
                    v_he[v1].compare_exchange_strong(
                        expected, base + 1, std::memory_order_relaxed);
                    expected = -1;
                    v_he[v2].compare_exchange_strong(
                        expected, base + 2, std::memory_order_relaxed);
                }
                };

            const size_t chunk = (nTri + nthreads - 1) / nthreads;
            for (unsigned i = 0; i < nthreads; ++i)
            {
                const size_t t0 = i * chunk;
                const size_t t1 = std::min(nTri, t0 + chunk);
                if (t0 >= t1) break;
                workers.emplace_back(fill, t0, t1);
            }
            for (auto& th : workers) th.join();
        }

        for (int i = 0; i < nV; ++i)
            mesh_vertices_[i].halfedge = v_he[i].load(std::memory_order_relaxed);

        std::vector<int> offset(nV + 1, 0);
        {
            MeshHalfEdge* hp = mesh_halfedges_.data();
            for (size_t i = 0; i < nHE; ++i)
                ++offset[hp[i].origin + 1];
            for (int v = 0; v < nV; ++v)
                offset[v + 1] += offset[v];
        }

        std::vector<int> adj(nHE);
        {
            MeshHalfEdge* hp = mesh_halfedges_.data();
            std::vector<int> cursor(offset.begin(), offset.begin() + nV);
            for (size_t i = 0; i < nHE; ++i)
            {
                const int u = hp[i].origin;
                adj[cursor[u]++] = static_cast<int>(i);
            }
        }

        {
            std::vector<std::thread> workers;
            workers.reserve(nthreads);
            ThreadJoiner joiner(workers);

            auto find_twins = [&](size_t i0, size_t i1) {
                MeshHalfEdge* hp = mesh_halfedges_.data();
                for (size_t i = i0; i < i1; ++i)
                {
                    if (hp[i].twin >= 0) continue;

                    const int u = hp[i].origin;
                    const int v = hp[hp[i].next].origin;

                    const int start = offset[v];
                    const int end = offset[v + 1];
                    for (int k = start; k < end; ++k)
                    {
                        const int h2 = adj[k];
                        if (hp[h2].origin != v) continue;
                        const int u2 = hp[hp[h2].next].origin;
                        if (u2 == u)
                        {
                            hp[i].twin = h2;
                            break;
                        }
                    }
                }
                };

            const size_t chunk = (nHE + nthreads - 1) / nthreads;
            for (unsigned i = 0; i < nthreads; ++i)
            {
                const size_t i0 = i * chunk;
                const size_t i1 = std::min(nHE, i0 + chunk);
                if (i0 >= i1) break;
                workers.emplace_back(find_twins, i0, i1);
            }
            for (auto& th : workers) th.join();
        }

        orig_to_local_ = std::move(orig_to_local);
        mesh_built_ = true;
    }

    int Triangulation::alloc_halfedge_()
    {
        if (!free_halfedges_.empty())
        {
            const int h = free_halfedges_.back();
            free_halfedges_.pop_back();
            return h;
        }
        const int h = static_cast<int>(mesh_halfedges_.size());
        mesh_halfedges_.push_back({});
        return h;
    }

    void Triangulation::free_halfedge_(int h)
    {
        mesh_halfedges_[h].origin = -1;
        mesh_halfedges_[h].twin = -1;
        mesh_halfedges_[h].next = -1;
        mesh_halfedges_[h].triangle_idx = -1;
        mesh_halfedges_[h].is_constrained = false;
        free_halfedges_.push_back(h);
    }

    int Triangulation::locate_(const qdt_vec2d& p, int& hull_he)
    {
        hull_he = -1;
        if (mesh_halfedges_.empty()) return -1;

        int h = last_face_he_;
        if (h < 0 || h >= static_cast<int>(mesh_halfedges_.size()) ||
            mesh_halfedges_[h].origin < 0)
        {
            h = 0;
            while (h < static_cast<int>(mesh_halfedges_.size()) &&
                mesh_halfedges_[h].origin < 0)
                ++h;
            if (h >= static_cast<int>(mesh_halfedges_.size())) return -1;
        }

        // A walk that is not stuck visits every triangle at most once. A
        // triangulation that is not Delaunay (constrained edges that were
        // flipped into place) can make the plain visibility walk run in
        // circles, so after this many steps the point is looked up by
        // scanning instead.
        const size_t max_steps = triangle_to_he_.size() + 16;
        for (size_t step = 0; step < max_steps; ++step)
        {
            const int h0 = h;
            const int h1 = mesh_halfedges_[h0].next;
            const int h2 = mesh_halfedges_[h1].next;
            if (mesh_halfedges_[h2].next != h0) break;

            const int v0 = mesh_halfedges_[h0].origin;
            const int v1 = mesh_halfedges_[h1].origin;
            const int v2 = mesh_halfedges_[h2].origin;

            const qdt_vec2d& A = mesh_vertices_[v0].p;
            const qdt_vec2d& B = mesh_vertices_[v1].p;
            const qdt_vec2d& C = mesh_vertices_[v2].p;

            const double o0 = robust::orient2d(A, B, p);
            const double o1 = robust::orient2d(B, C, p);
            const double o2 = robust::orient2d(C, A, p);

            if (o0 >= 0 && o1 >= 0 && o2 >= 0)
            {
                last_face_he_ = h0;
                return h0;
            }

            int cross;
            if (o0 < 0)      cross = h0;
            else if (o1 < 0) cross = h1;
            else             cross = h2;

            const int tw = mesh_halfedges_[cross].twin;
            if (tw < 0)
            {
                hull_he = cross;
                return -1;
            }

            h = tw;
        }

        return locate_scan_(p, hull_he);
    }

    // Finds the triangle that contains p (boundary included) by testing all of
    // them. Returns a halfedge of that triangle. If p is outside the convex
    // hull, returns -1 and sets hull_he to a hull halfedge that p can see.
    int Triangulation::locate_scan_(const qdt_vec2d& p, int& hull_he)
    {
        hull_he = -1;
        for (size_t t = 0; t < triangle_to_he_.size(); ++t)
        {
            const int h0 = triangle_to_he_[t];
            const int h1 = mesh_halfedges_[h0].next;
            const int h2 = mesh_halfedges_[h1].next;
            const qdt_vec2d& A = mesh_vertices_[mesh_halfedges_[h0].origin].p;
            const qdt_vec2d& B = mesh_vertices_[mesh_halfedges_[h1].origin].p;
            const qdt_vec2d& C = mesh_vertices_[mesh_halfedges_[h2].origin].p;
            if (robust::orient2d(A, B, p) >= 0 && robust::orient2d(B, C, p) >= 0 &&
                robust::orient2d(C, A, p) >= 0)
            {
                last_face_he_ = h0;
                return h0;
            }
        }
        for (size_t h = 0; h < mesh_halfedges_.size(); ++h)
        {
            const MeshHalfEdge& e = mesh_halfedges_[h];
            if (e.origin < 0 || e.twin >= 0) continue;
            if (robust::orient2d(mesh_vertices_[e.origin].p,
                mesh_vertices_[mesh_halfedges_[e.next].origin].p, p) < 0)
            {
                hull_he = static_cast<int>(h);
                return -1;
            }
        }
        return -1;
    }

    void Triangulation::insert_by_rebuild_(const qdt_vec2d& p)
    {
        input_points.push_back(p);
        Triangulation fresh = delaunay_triangulate(
            input_points.data(), input_points.size());
        *this = std::move(fresh);
    }

    void Triangulation::insert(const qdt_vec2d& p)
    {
        ensure_mesh_();
        insert_after_mesh_(p);
    }

    void Triangulation::insert_parallel(const qdt_vec2d& p)
    {
        ensure_mesh_parallel_();
        insert_after_mesh_(p);
    }

    void Triangulation::insert_after_mesh_(const qdt_vec2d& p)
    {
        if (mesh_halfedges_.empty())
        {
            insert_by_rebuild_(p);
            return;
        }

        int hull_he = -1;
        const int h_start = locate_(p, hull_he);
        if (h_start < 0 && hull_he < 0)
        {
            insert_by_rebuild_(p);
            return;
        }
        const bool outside = (h_start < 0);

        // Duplicate point. locate_ returns a triangle that contains p, border
        // included, tested with exact predicates. A vertex of the triangulation
        // is never inside a triangle or inside one of its edges, so if p has the
        // coordinates of an existing vertex, that vertex is a corner of the
        // triangle found. Looking at the three corners is therefore enough; it
        // replaces a scan over all vertices. A duplicate cannot lie outside the
        // convex hull, so the outside case needs no check.
        if (!outside)
        {
            int e = h_start;
            do
            {
                const MeshVertex& v = mesh_vertices_[mesh_halfedges_[e].origin];
                if (v.p.x == p.x && v.p.y == p.y)
                {
                    input_points.push_back(p);
                    rep.push_back(v.orig_id);
                    return;
                }
                e = mesh_halfedges_[e].next;
            } while (e != h_start);
        }

        const size_t new_id = input_points.size();
        input_points.push_back(p);
        rep.push_back(new_id);
        ++unique_points;

        const int new_v = static_cast<int>(mesh_vertices_.size());
        {
            MeshVertex v;
            v.p = p;
            v.orig_id = new_id;
            v.halfedge = -1;
            mesh_vertices_.push_back(v);
        }
        // Keep the point -> mesh vertex table in step (see orig_to_local_).
        if (orig_to_local_.size() <= new_id)
            orig_to_local_.resize(new_id + 1, -1);
        orig_to_local_[new_id] = new_v;

        auto he_dest = [&](int e) {
            return mesh_halfedges_[mesh_halfedges_[e].next].origin;
            };

        std::vector<int> visible_hull;
        if (outside)
        {
            auto is_visible = [&](int e) {
                const qdt_vec2d& A = mesh_vertices_[mesh_halfedges_[e].origin].p;
                const qdt_vec2d& B = mesh_vertices_[he_dest(e)].p;
                return robust::orient2d(A, B, p) < 0;
                };
            auto next_hull = [&](int e) {
                int x = mesh_halfedges_[e].next;
                while (mesh_halfedges_[x].twin >= 0)
                    x = mesh_halfedges_[mesh_halfedges_[x].twin].next;
                return x;
                };
            auto prev_hull = [&](int e) {
                int x = mesh_halfedges_[mesh_halfedges_[e].next].next;
                while (mesh_halfedges_[x].twin >= 0)
                {
                    const int tw = mesh_halfedges_[x].twin;
                    x = mesh_halfedges_[mesh_halfedges_[tw].next].next;
                }
                return x;
                };

            visible_hull.push_back(hull_he);
            for (int e = next_hull(hull_he); e != hull_he && is_visible(e); e = next_hull(e))
                visible_hull.push_back(e);
            for (int e = prev_hull(hull_he); e != hull_he && is_visible(e); e = prev_hull(e))
                visible_hull.push_back(e);
        }

        // The cavity is only the triangle that contains p, plus the triangle on
        // the other side of the edge p lies on (if any), so that edge is split
        // and no degenerate triangle appears. Growing the cavity with
        // circumcircle tests alone would also swallow constrained edges inside
        // the grown region and delete them. Delaunay is restored afterwards by
        // edge flips (see the end of this function), which never flip a
        // constrained edge, so constraints survive.
        // If p is outside the hull nothing is removed: p is joined to the
        // visible hull edges and the flips take care of the rest.
        std::vector<int>  cavity_hes;
        // Ids of the triangles in the cavity: at most two (the one containing p
        // and, if p lies on an edge, the one across it), so a linear search is
        // cheaper than a flag per triangle of the whole triangulation.
        std::vector<int>  cavity_tris;
        auto in_cavity = [&](int t) {
            return std::find(cavity_tris.begin(), cavity_tris.end(), t) != cavity_tris.end();
            };
        int split_a = -1, split_b = -1; // endpoints of the constrained edge p splits

        if (!outside)
        {
            auto add_to_cavity = [&](int h) {
                const int t = mesh_halfedges_[h].triangle_idx;
                if (t >= 0 && !in_cavity(t))
                {
                    cavity_tris.push_back(t);
                    cavity_hes.push_back(h);
                }
                };

            add_to_cavity(h_start);
            int e = h_start;
            do
            {
                const int e_next = mesh_halfedges_[e].next;
                if (robust::orient2d(mesh_vertices_[mesh_halfedges_[e].origin].p,
                    mesh_vertices_[mesh_halfedges_[e_next].origin].p, p) == 0)
                {
                    // p lies on this edge, so the triangle across it is split too.
                    const int tw = mesh_halfedges_[e].twin;
                    if (tw >= 0) add_to_cavity(tw);
                    if (mesh_halfedges_[e].is_constrained)
                    {
                        split_a = mesh_halfedges_[e].origin;
                        split_b = mesh_halfedges_[e_next].origin;
                    }
                }
                e = e_next;
            } while (e != h_start);
        }

        struct BoundaryInfo { int a, b, external_twin; bool is_constrained; };
        std::vector<BoundaryInfo> boundary_info;

        for (int h : cavity_hes)
        {
            int e = h;
            do
            {
                const int tw = mesh_halfedges_[e].twin;
                int tw_tri = (tw >= 0) ? mesh_halfedges_[tw].triangle_idx : -1;
                bool tw_in_cavity = (tw_tri >= 0 && in_cavity(tw_tri));

                if (!tw_in_cavity)
                {
                    const int a = mesh_halfedges_[e].origin;
                    const int b = mesh_halfedges_[mesh_halfedges_[e].next].origin;
                    if (tw >= 0 ||
                        robust::orient2d(mesh_vertices_[a].p, mesh_vertices_[b].p, p) > 0)
                        boundary_info.push_back({ a, b, tw, mesh_halfedges_[e].is_constrained });
                }

                e = mesh_halfedges_[e].next;
            } while (e != h);
        }

        for (int e : visible_hull)
        {
            const int tri_idx = mesh_halfedges_[e].triangle_idx;
            if (tri_idx >= 0 && in_cavity(tri_idx))
                continue;
            const int e1 = mesh_halfedges_[e].next;
            const int a = mesh_halfedges_[e].origin;
            const int b = mesh_halfedges_[e1].origin;
            boundary_info.push_back({ b, a, e, mesh_halfedges_[e].is_constrained });
        }

        // Highest id first: remove_triangle_by_index_ moves the last triangle
        // into the freed slot, which must not invalidate an id still to come.
        std::sort(cavity_tris.rbegin(), cavity_tris.rend());
        for (int idx : cavity_tris)
        {
            remove_triangle_by_index_(idx);
        }

        for (int h : cavity_hes)
        {
            int e = h;
            do
            {
                const int nxt = mesh_halfedges_[e].next;
                free_halfedge_(e);
                e = nxt;
            } while (e != h);
        }

        struct NewTri { int h0, h1, h2; };
        std::vector<NewTri> new_tris;

        for (const auto& bi : boundary_info)
        {
            const int h0 = alloc_halfedge_();
            const int h1 = alloc_halfedge_();
            const int h2 = alloc_halfedge_();

            mesh_halfedges_[h0] = { bi.a, bi.external_twin, h1, -1, bi.is_constrained };
            mesh_halfedges_[h1] = { bi.b, -1, h2, -1, false };
            mesh_halfedges_[h2] = { new_v, -1, h0, -1, false };

            if (bi.external_twin >= 0)
            {
                mesh_halfedges_[bi.external_twin].twin = h0;
                mesh_halfedges_[bi.external_twin].is_constrained = bi.is_constrained;
            }

            size_t new_tri_idx = add_triangle_to_output_(
                mesh_vertices_[bi.a].orig_id,
                mesh_vertices_[bi.b].orig_id,
                new_id,
                h0
            );

            mesh_halfedges_[h0].triangle_idx = static_cast<int>(new_tri_idx);
            mesh_halfedges_[h1].triangle_idx = static_cast<int>(new_tri_idx);
            mesh_halfedges_[h2].triangle_idx = static_cast<int>(new_tri_idx);

            new_tris.push_back({ h0, h1, h2 });
        }

        if (new_tris.empty())
        {
            Triangulation fresh = delaunay_triangulate(
                input_points.data(), input_points.size());
            *this = std::move(fresh);
            return;
        }

        std::unordered_map<int, int> v_to_newv;
        std::unordered_map<int, int> newv_to_v;

        for (const auto& nt : new_tris)
        {
            const int a = mesh_halfedges_[nt.h0].origin;
            const int b = mesh_halfedges_[nt.h1].origin;
            v_to_newv[b] = nt.h1;
            newv_to_v[a] = nt.h2;
        }

        for (const auto& kv : v_to_newv)
        {
            const int v = kv.first;
            const int h1 = kv.second;
            auto it = newv_to_v.find(v);
            if (it != newv_to_v.end())
            {
                const int h2 = it->second;
                mesh_halfedges_[h1].twin = h2;
                mesh_halfedges_[h2].twin = h1;
            }
        }

        mesh_vertices_[new_v].halfedge = new_tris[0].h2;
        last_face_he_ = new_tris[0].h0;

        // The removed triangles' halfedges were freed and their slots reused
        // above, so the stored halfedge of every cavity vertex may be stale.
        // h0 starts at boundary_info[i].a and h1 at boundary_info[i].b.
        for (size_t i = 0; i < new_tris.size(); ++i)
        {
            mesh_vertices_[boundary_info[i].a].halfedge = new_tris[i].h0;
            mesh_vertices_[boundary_info[i].b].halfedge = new_tris[i].h1;
        }

        // p split a constrained edge (split_a, split_b): both halves stay
        // constrained. new_tris[i] was built from boundary_info[i]; its h1
        // runs b -> p and its h2 runs p -> a.
        if (split_a >= 0)
        {
            auto mark_constrained = [&](int h) {
                mesh_halfedges_[h].is_constrained = true;
                const int tw = mesh_halfedges_[h].twin;
                if (tw >= 0) mesh_halfedges_[tw].is_constrained = true;
                };
            for (size_t i = 0; i < new_tris.size(); ++i)
            {
                if (boundary_info[i].b == split_a || boundary_info[i].b == split_b)
                    mark_constrained(new_tris[i].h1);
                if (boundary_info[i].a == split_a || boundary_info[i].a == split_b)
                    mark_constrained(new_tris[i].h2);
            }
        }

        // Restore the (constrained) Delaunay property around p with edge
        // flips. h0 of every new triangle is the edge opposite p; flipping it
        // puts a diagonal through p, and the two edges that then face p are
        // checked in turn. Constrained edges are never flipped.
        std::vector<int> pending;
        for (const auto& nt : new_tris)
            pending.push_back(nt.h0);
        while (!pending.empty())
        {
            const int h = pending.back();
            pending.pop_back();

            const int tw = mesh_halfedges_[h].twin;
            if (tw < 0 || mesh_halfedges_[h].is_constrained)
                continue;

            const int h_next = mesh_halfedges_[h].next;
            const int h_prev = mesh_halfedges_[h_next].next;
            if (mesh_halfedges_[h_prev].origin != new_v)
                continue; // h's triangle no longer has p as its apex

            const int tw_next = mesh_halfedges_[tw].next;
            const int tw_prev = mesh_halfedges_[tw_next].next;

            const int a = mesh_halfedges_[h].origin;
            const int b = mesh_halfedges_[h_next].origin;
            const int d = mesh_halfedges_[tw_prev].origin;

            if (robust::incircle(mesh_vertices_[a].p, mesh_vertices_[b].p,
                mesh_vertices_[new_v].p, mesh_vertices_[d].p) <= 0)
                continue;

            if (flip_edge_(h))
            {
                pending.push_back(tw_next); // a -> d
                pending.push_back(tw_prev); // d -> b
            }
        }
    }

    int Triangulation::orig_to_local_vtx_(size_t orig_id) const
    {
        // Constant time: the table is maintained by ensure_mesh_*() and
        // insert_after_mesh_(). Ids past its end have no mesh vertex.
        if (orig_id >= orig_to_local_.size())
            return -1;
        return orig_to_local_[orig_id];
    }

    // Calls f(h) for every halfedge h that starts at vertex x, until f returns
    // true. Works for interior and for convex hull vertices (whose fan is open
    // on the outer side), independent of which halfedge is stored for x.
    template <class F>
    void Triangulation::for_each_out_(int x, F&& f) const
    {
        const int h0 = mesh_vertices_[x].halfedge;
        if (h0 < 0) return;
        const size_t limit = mesh_halfedges_.size();
        size_t steps = 0;

        // Counter-clockwise from h0: the next halfedge out of x is the twin of
        // the previous halfedge of the same triangle.
        int h = h0;
        for (;;)
        {
            if (f(h)) return;
            if (++steps > limit) return;
            const int prev = mesh_halfedges_[mesh_halfedges_[h].next].next;
            const int tw = mesh_halfedges_[prev].twin;
            if (tw < 0) break;      // reached the hull boundary
            if (tw == h0) return;   // went all the way around
            h = tw;
        }

        // Clockwise from h0 (only reached for a hull vertex).
        h = h0;
        for (;;)
        {
            const int tw = mesh_halfedges_[h].twin;
            if (tw < 0) return;
            h = mesh_halfedges_[tw].next;
            if (h == h0) return;
            if (f(h)) return;
            if (++steps > limit) return;
        }
    }

    // The halfedge that goes from u to v, or -1. A convex hull edge has a
    // halfedge in one direction only, so callers that accept either direction
    // look for (u, v) and then for (v, u).
    int Triangulation::find_halfedge_(int u_local, int v_local) const
    {
        int found = -1;
        for_each_out_(u_local, [&](int h) -> bool {
            if (mesh_halfedges_[mesh_halfedges_[h].next].origin == v_local)
            {
                found = h;
                return true;
            }
            return false;
            });
        return found;
    }

    bool Triangulation::flip_edge_(int h)
    {
        int tw = mesh_halfedges_[h].twin;
        if (tw < 0) return false;
        if (mesh_halfedges_[h].is_constrained) return false;

        int h_next = mesh_halfedges_[h].next;
        int h_prev = mesh_halfedges_[h_next].next;

        int tw_next = mesh_halfedges_[tw].next;
        int tw_prev = mesh_halfedges_[tw_next].next;

        // Before the flip (both triangles counter-clockwise):
        //   triangle 1:  h : a->b,  h_next : b->c,  h_prev : c->a
        //   triangle 2:  tw: b->a,  tw_next: a->d,  tw_prev: d->b
        // so c is the apex of triangle 1 and d is the apex of triangle 2.
        // (tw originates at b and tw_next originates at a, so neither can be
        // used for c or d.)
        const int a = mesh_halfedges_[h].origin;
        const int b = mesh_halfedges_[h_next].origin;
        const int c = mesh_halfedges_[h_prev].origin;
        const int d = mesh_halfedges_[tw_prev].origin;

        // After the flip the diagonal is c-d and the triangles are (d,c,a) and
        // (c,d,b). The flip is only valid when both stay strictly
        // counter-clockwise, i.e. the quadrilateral a,d,b,c is strictly convex.
        if (robust::orient2d(mesh_vertices_[d].p, mesh_vertices_[c].p, mesh_vertices_[a].p) <= 0 ||
            robust::orient2d(mesh_vertices_[c].p, mesh_vertices_[d].p, mesh_vertices_[b].p) <= 0)
        {
            return false;
        }

        // triangle 1 becomes  h : d->c,  h_prev : c->a,  tw_next : a->d
        mesh_halfedges_[h].origin = d;
        mesh_halfedges_[h].next = h_prev;
        mesh_halfedges_[h_prev].next = tw_next;
        mesh_halfedges_[tw_next].next = h;

        // triangle 2 becomes  tw : c->d,  tw_prev : d->b,  h_next : b->c
        mesh_halfedges_[tw].origin = c;
        mesh_halfedges_[tw].next = tw_prev;
        mesh_halfedges_[tw_prev].next = h_next;
        mesh_halfedges_[h_next].next = tw;

        const int tri1 = mesh_halfedges_[h].triangle_idx;
        const int tri2 = mesh_halfedges_[tw].triangle_idx;

        // h and h_prev stay in triangle 1, tw and tw_prev stay in triangle 2;
        // tw_next moves to triangle 1 and h_next moves to triangle 2.
        mesh_halfedges_[tw_next].triangle_idx = tri1;
        mesh_halfedges_[h_next].triangle_idx = tri2;

        triangles[tri1] = { mesh_vertices_[d].orig_id, mesh_vertices_[c].orig_id, mesh_vertices_[a].orig_id };
        triangles[tri2] = { mesh_vertices_[c].orig_id, mesh_vertices_[d].orig_id, mesh_vertices_[b].orig_id };

        // remove_triangle_by_index_ walks the cycle starting at
        // triangle_to_he_[i], so the representative must lie in triangle i.
        triangle_to_he_[tri1] = h;
        triangle_to_he_[tri2] = tw;

        // A vertex's stored halfedge must start at that vertex.
        mesh_vertices_[a].halfedge = tw_next;
        mesh_vertices_[b].halfedge = h_next;
        mesh_vertices_[c].halfedge = h_prev;
        mesh_vertices_[d].halfedge = h;

        return true;
    }

    bool Triangulation::segment_intersects_(const qdt_vec2d& a, const qdt_vec2d& b, const qdt_vec2d& c, const qdt_vec2d& d) const
    {
        double o1 = robust::orient2d(a, b, c);
        double o2 = robust::orient2d(a, b, d);
        double o3 = robust::orient2d(c, d, a);
        double o4 = robust::orient2d(c, d, b);

        return ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
            ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0));
    }

    void Triangulation::mark_constrained_(int h)
    {
        mesh_halfedges_[h].is_constrained = true;
        const int tw = mesh_halfedges_[h].twin;
        if (tw >= 0) mesh_halfedges_[tw].is_constrained = true;
    }

    // Do the vectors a->b and a->c point the same way? Only meaningful when
    // they are known to be collinear. Compares the signs of the coordinate
    // differences, so it cannot underflow.
    static bool same_direction_(const qdt_vec2d& a, const qdt_vec2d& b, const qdt_vec2d& c)
    {
        auto sg = [](double d) { return d > 0 ? 1 : (d < 0 ? -1 : 0); };
        const int sx = sg(b.x - a.x) * sg(c.x - a.x);
        const int sy = sg(b.y - a.y) * sg(c.y - a.y);
        return sx >= 0 && sy >= 0 && (sx > 0 || sy > 0);
    }

    // Walks from u to v and lists, in order, the edges the open segment uv
    // crosses. Each is returned as the halfedge that runs from the vertex on
    // the right of u->v to the vertex on its left. Nothing is modified.
    //
    // Returns QDT_OK, or (and the list is meaningless) the reason why the
    // segment cannot become an edge as it is: QDT_ERR_CROSSES_CONSTRAINT, or
    // QDT_ERR_THROUGH_POINT when it passes exactly through another vertex.
    // QDT_ERR_INTERNAL is returned for situations that cannot occur when both
    // end points are vertices of the triangulation (the walk leaves the convex
    // hull, or does not end).
    qdt_status Triangulation::trace_segment_(int u, int v, std::vector<int>& crossed) const
    {
        auto he_next = [&](int h) { return mesh_halfedges_[h].next; };
        auto org = [&](int h) { return mesh_halfedges_[h].origin; };
        auto pt = [&](int x) -> const qdt_vec2d& { return mesh_vertices_[x].p; };
        const qdt_vec2d& U = pt(u);
        const qdt_vec2d& V = pt(v);

        // Among the triangles around u, find the one the segment leaves u
        // through: v is strictly between its edges u->y and u->z. If v lies
        // exactly on one of those edges, that edge's other end is a vertex on
        // the segment.
        int start = -1; // y -> z, the edge opposite u in that triangle
        bool blocked = false;
        for_each_out_(u, [&](int h) -> bool {
            const int hn = he_next(h);
            const int y = org(hn);
            const int z = org(he_next(hn));
            const double oy = robust::orient2d(U, pt(y), V);
            const double oz = robust::orient2d(U, pt(z), V);
            if (oy > 0 && oz < 0)
            {
                start = hn;
                return true;
            }
            if ((oy == 0 && same_direction_(U, pt(y), V)) ||
                (oz == 0 && same_direction_(U, pt(z), V)))
            {
                blocked = true;
                return true;
            }
            return false;
            });
        if (blocked) return QDT_ERR_THROUGH_POINT;
        if (start < 0) return QDT_ERR_INTERNAL; // v is a vertex, so some triangle must contain the direction

        // In the first triangle y lies on the right of u->v and z on its left.
        // Every crossed edge keeps that orientation (right vertex -> left vertex).
        int c = start;
        for (size_t guard = 0; guard <= triangles.size(); ++guard)
        {
            if (mesh_halfedges_[c].is_constrained) return QDT_ERR_CROSSES_CONSTRAINT;
            crossed.push_back(c);

            const int tw = mesh_halfedges_[c].twin; // z -> y in the next triangle
            if (tw < 0) return QDT_ERR_INTERNAL;    // would leave the convex hull
            const int a = he_next(tw);              // y -> w
            const int b = he_next(a);               // w -> z
            const int w = org(b);
            if (w == v) return QDT_OK;

            const double o = robust::orient2d(U, V, pt(w));
            if (o == 0) return QDT_ERR_THROUGH_POINT; // w lies on the segment
            c = (o > 0) ? a : b;                    // w on the left: leave via y->w
        }
        return QDT_ERR_INTERNAL;
    }

    // Lawson's flip algorithm: every edge on the stack that is not
    // constrained and not locally Delaunay is flipped, and the four edges
    // around the flipped one are examined in turn.
    void Triangulation::legalize_(std::vector<int>& stack)
    {
        while (!stack.empty())
        {
            const int h = stack.back();
            stack.pop_back();

            if (mesh_halfedges_[h].origin < 0) continue;
            const int tw = mesh_halfedges_[h].twin;
            if (tw < 0 || mesh_halfedges_[h].is_constrained) continue;

            const int h_next = mesh_halfedges_[h].next;
            const int h_prev = mesh_halfedges_[h_next].next;
            const int tw_next = mesh_halfedges_[tw].next;
            const int tw_prev = mesh_halfedges_[tw_next].next;

            const int a = mesh_halfedges_[h].origin;
            const int b = mesh_halfedges_[h_next].origin;
            const int c = mesh_halfedges_[h_prev].origin;
            const int d = mesh_halfedges_[tw_prev].origin;

            if (robust::incircle(mesh_vertices_[a].p, mesh_vertices_[b].p,
                mesh_vertices_[c].p, mesh_vertices_[d].p) <= 0)
                continue;
            if (!flip_edge_(h)) continue;

            stack.push_back(h_prev);
            stack.push_back(h_next);
            stack.push_back(tw_next);
            stack.push_back(tw_prev);
        }
    }

    // Inserts the constraint between two mesh vertices (local indices).
    //  - If the edge exists it just becomes constrained.
    //  - Otherwise the edges the segment crosses are flipped away (Sloan 1993:
    //    a crossed edge whose quadrilateral is strictly convex is flipped; if
    //    the new diagonal still crosses the segment it goes back in the queue),
    //    which always ends with the edge u-v in place, and then the
    //    constrained Delaunay property is restored with Lawson flips.
    //  - If the segment crosses a constrained edge or passes through another
    //    vertex it cannot be an edge: nothing is changed and the reason
    //    (QDT_ERR_CROSSES_CONSTRAINT or QDT_ERR_THROUGH_POINT) is returned.
    //  - QDT_ERR_INTERNAL means the flipping did not end with the edge in
    //    place. That cannot happen with exact predicates; the triangulation may
    //    have been changed, so the caller resets the context.
    qdt_status Triangulation::insert_constraint_local_(int u, int v)
    {
        int e = find_halfedge_(u, v);
        if (e < 0) e = find_halfedge_(v, u);
        if (e >= 0)
        {
            mark_constrained_(e);
            return QDT_OK;
        }

        std::vector<int> crossed;
        const qdt_status traced = trace_segment_(u, v, crossed);
        if (traced != QDT_OK) return traced;

        const qdt_vec2d& U = mesh_vertices_[u].p;
        const qdt_vec2d& V = mesh_vertices_[v].p;

        std::deque<int> queue(crossed.begin(), crossed.end());
        std::vector<int> new_edges;
        size_t stalled = 0;
        while (!queue.empty())
        {
            const int h = queue.front();
            queue.pop_front();
            if (flip_edge_(h))
            {
                stalled = 0;
                const int c = mesh_halfedges_[h].origin;
                const int d = mesh_halfedges_[mesh_halfedges_[h].next].origin;
                if (segment_intersects_(U, V, mesh_vertices_[c].p, mesh_vertices_[d].p))
                    queue.push_back(h);
                else
                    new_edges.push_back(h);
            }
            else
            {
                queue.push_back(h);
                // A whole round without a single flip. With no vertex on the
                // segment this cannot happen; stop instead of looping.
                if (++stalled >= queue.size()) break;
            }
        }

        bool ok = false;
        if (queue.empty())
        {
            e = -1;
            for (int h : new_edges)
            {
                const int a = mesh_halfedges_[h].origin;
                const int b = mesh_halfedges_[mesh_halfedges_[h].next].origin;
                if ((a == u && b == v) || (a == v && b == u))
                {
                    e = h;
                    break;
                }
            }
            if (e < 0) e = find_halfedge_(u, v);
            if (e < 0) e = find_halfedge_(v, u);
            if (e >= 0)
            {
                mark_constrained_(e);
                ok = true;
            }
        }

        // Restore the Delaunay property around everything that changed. The
        // constrained edge u-v is already marked, so it is never flipped.
        std::vector<int> stack;
        auto push_around = [&](int h) {
            if (mesh_halfedges_[h].origin < 0) return;
            const int n = mesh_halfedges_[h].next;
            stack.push_back(h);
            stack.push_back(n);
            stack.push_back(mesh_halfedges_[n].next);
            const int tw = mesh_halfedges_[h].twin;
            if (tw >= 0)
            {
                const int tn = mesh_halfedges_[tw].next;
                stack.push_back(tn);
                stack.push_back(mesh_halfedges_[tn].next);
            }
            };
        for (int h : new_edges) push_around(h);
        for (int h : queue) push_around(h);
        legalize_(stack);
        return ok ? QDT_OK : QDT_ERR_INTERNAL;
    }

    qdt_status Triangulation::insert_constraint(size_t orig_v0, size_t orig_v1)
    {
        // Indices come from the caller: reject anything that is not a point of
        // this context (this also covers a context that was never filled).
        if (orig_v0 >= rep.size() || orig_v1 >= rep.size()) return QDT_ERR_INVALID_INDEX;
        orig_v0 = rep[orig_v0];
        orig_v1 = rep[orig_v1];
        if (orig_v0 == orig_v1) return QDT_ERR_SAME_POINT;

        ensure_mesh_();
        const int u = orig_to_local_vtx_(orig_v0);
        const int v = orig_to_local_vtx_(orig_v1);
        // Every distinct point is a vertex of the mesh as soon as there are
        // triangles, so a missing vertex means there are none.
        if (u < 0 || v < 0) return QDT_ERR_NO_TRIANGLES;
        return insert_constraint_local_(u, v);
    }

    qdt_status Triangulation::insert_constraints_batch(const qdt_edge* edges, size_t num_edges, qdt_status* results)
    {
        qdt_status first_error = QDT_OK;
        for (size_t i = 0; i < num_edges; ++i)
        {
            const qdt_status st = insert_constraint(edges[i].v0, edges[i].v1);
            if (st == QDT_ERR_INTERNAL)
            {
                // The context is no longer trustworthy: stop, and let the
                // caller reset it. No entry counts as accepted.
                if (results)
                {
                    for (size_t k = 0; k < num_edges; ++k)
                        results[k] = QDT_ERR_INTERNAL;
                }
                return QDT_ERR_INTERNAL;
            }
            if (results) results[i] = st;
            if (st != QDT_OK && first_error == QDT_OK) first_error = st;
        }
        return first_error;
    }

    Triangulation delaunay_triangulate(const Point* pts, size_t sz)
    {
        Triangulation res;
        std::vector<detail::Vtx> V;
        detail::load_unique_points(pts, sz, res, V);

        if (V.size() >= 3)
        {
            detail::Builder<detail::SerialPolicy> b(V);
            b.run();
            b.extract(res.triangles);
        }

        res.input_points.resize(sz);
        for (size_t i = 0; i < sz; ++i)
            res.input_points[i] = pts[i];

        return res;
    }
}