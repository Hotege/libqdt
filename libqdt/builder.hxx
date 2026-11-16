/*
 * Copyright (c) 2026 QDT contributors
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the QDT library.
 */
 // The divide-and-conquer Delaunay builder, shared by the serial implementation
 // (GuibasStolfiDwyer_Wrapper.cpp) and the parallel one
 // (GuibasStolfiDwyer_Threads.cpp). Internal header: it is not part of the
 // public interface and is not installed.
 //
 // The algorithm is Guibas and Stolfi's divide and conquer on the quad-edge
 // structure, with Dwyer's alternating vertical and horizontal cuts. It is
 // written once, as a template over a Policy that decides two things the two
 // implementations do differently:
 //
 //   * where the quad-edge records live (a growing vector with a free list, or
 //     blocks of a pre-allocated array handed out per thread), and
 //   * whether the two halves of a build, or the two scans of a merge step, run
 //     one after the other or on different threads.
 //
 // A Policy must be default-constructible and provide:
 //
 //   void init(std::vector<QE>& qe, size_t nverts);
 //       Called once by the constructor, to size the edge storage.
 //   int  alloc(std::vector<QE>& qe);
 //       Returns the index of an unused quad-edge record (a free-list entry
 //       first, otherwise a new one). May throw.
 //   void release(std::vector<QE>& qe, int k);
 //       Gives record k back; it must be reusable by a later alloc().
 //   int  edge_count(const std::vector<QE>& qe) const;
 //       The number of records extract() has to look at: every record that
 //       alloc() has ever returned has an index below it.
 //   bool split_build(int depth, int n) const;
 //       Whether the two halves of a build of n points at this depth may run
 //       concurrently.
 //   bool split_merge(int depth) const;
 //       Whether the two scans of a merge step at this depth may run
 //       concurrently.
 //   auto fork_join(bool concurrent, A&& first, B&& second)
 //       -> std::pair<decltype(first()), decltype(second())>;
 //       Calls first() and second() and returns both results. With
 //       concurrent == false, first() runs to completion before second()
 //       starts, on the calling thread. With true, first() may run on another
 //       thread while second() runs on the calling thread; the call returns
 //       when both are done. first must capture by value: if second() throws,
 //       the call unwinds while first() may still be running.
#if !defined(QDT_BUILDER_HXX_)
#define QDT_BUILDER_HXX_

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

#include "qdt.h"
#include "robust.hxx"
#include "wrapper.hxx"

namespace wrapper
{
    namespace detail
    {
        struct Vtx : Point
        {
            size_t id;
        };

        inline bool lessF(const Point& a, const Point& b, int frame)
        {
            if (frame == 0)
                return a.x < b.x || (a.x == b.x && a.y < b.y);
            return a.y < b.y || (a.y == b.y && a.x > b.x);
        }

        inline int rot(int e) { return (e & ~3) | ((e + 1) & 3); }
        inline int invrot(int e) { return (e & ~3) | ((e + 3) & 3); }
        inline int sym(int e) { return e ^ 2; }

        struct QE
        {
            int next[4];
            int vtx[2];
        };

        template <class Policy>
        class Builder
        {
        public:
            explicit Builder(std::vector<Vtx>& v) : V(v) { policy_.init(qe, V.size()); }

            Policy& policy() { return policy_; }

            void run() { build(0, static_cast<int>(V.size()), 0, 0, true, 0); }

            void extract(std::vector<std::array<size_t, 3>>& out) const
            {
                out.reserve(2 * V.size());
                const int sz = policy_.edge_count(qe);
                for (int k = 0; k < sz; ++k)
                {
                    if (qe[k].vtx[0] < 0) continue;
                    for (int r = 0; r < 4; r += 2)
                    {
                        const int e = 4 * k + r;
                        const int f = lnext(e), g = lnext(f);
                        if (lnext(g) != e) continue;
                        if (e > f || e > g) continue;
                        const int a = org(e), b = org(f), c = org(g);
                        if (robust::orient2d(V[a], V[b], V[c]) > 0)
                            out.push_back({ V[a].id, V[b].id, V[c].id });
                    }
                }
            }

        private:
            std::vector<Vtx>& V;
            std::vector<QE> qe;
            Policy policy_;

            int onext(int e) const { return qe[e >> 2].next[e & 3]; }
            int org(int e) const { return qe[e >> 2].vtx[(e >> 1) & 1]; }
            int dest(int e) const { return org(sym(e)); }
            void setOrg(int e, int v) { qe[e >> 2].vtx[(e >> 1) & 1] = v; }
            int oprev(int e) const { return rot(onext(rot(e))); }
            int lnext(int e) const { return rot(onext(invrot(e))); }
            int rprev(int e) const { return onext(sym(e)); }

            int makeEdge()
            {
                const int k = policy_.alloc(qe);
                QE& q = qe[k];
                const int e = 4 * k;
                q.next[0] = e;
                q.next[1] = e + 3;
                q.next[2] = e + 2;
                q.next[3] = e + 1;
                q.vtx[0] = q.vtx[1] = -1;
                return e;
            }

            void splice(int a, int b)
            {
                const int alpha = rot(onext(a)), beta = rot(onext(b));
                std::swap(qe[a >> 2].next[a & 3], qe[b >> 2].next[b & 3]);
                std::swap(qe[alpha >> 2].next[alpha & 3], qe[beta >> 2].next[beta & 3]);
            }

            int connect(int a, int b)
            {
                const int e = makeEdge();
                setOrg(e, dest(a));
                setOrg(sym(e), org(b));
                splice(e, lnext(a));
                splice(sym(e), b);
                return e;
            }

            void deleteEdge(int e)
            {
                splice(e, oprev(e));
                splice(sym(e), oprev(sym(e)));
                policy_.release(qe, e >> 2);
            }

            bool ccw(int a, int b, int c) const { return robust::orient2d(V[a], V[b], V[c]) > 0; }
            bool rightOf(int x, int e) const { return ccw(x, dest(e), org(e)); }
            bool leftOf(int x, int e) const { return ccw(x, org(e), dest(e)); }
            bool inCircle(int a, int b, int c, int d) const
            {
                if (a == d || b == d || c == d || a == b || a == c || b == c) return false;
                return robust::incircle(V[a], V[b], V[c], V[d]) > 0;
            }

            std::pair<int, int> leaf(int lo, int hi, int frame)
            {
                std::sort(V.begin() + lo, V.begin() + hi,
                    [frame](const Vtx& a, const Vtx& b) { return lessF(a, b, frame); });
                const int n = hi - lo;
                const int s1 = lo, s2 = lo + 1;
                if (n == 2)
                {
                    const int a = makeEdge();
                    setOrg(a, s1);
                    setOrg(sym(a), s2);
                    return { a, sym(a) };
                }
                const int s3 = lo + 2;
                const int a = makeEdge(), b = makeEdge();
                splice(sym(a), b);
                setOrg(a, s1);
                setOrg(sym(a), s2);
                setOrg(b, s2);
                setOrg(sym(b), s3);
                const double o = robust::orient2d(V[s1], V[s2], V[s3]);
                if (o > 0) { connect(b, a); return { a, sym(b) }; }
                if (o < 0) { const int c = connect(b, a); return { sym(c), c }; }
                return { a, sym(b) };
            }

            // The two scans of a merge step. They only touch the edges on their
            // own side of the base edge, so they may run concurrently.
            std::pair<int, bool> merge_scan_left(int basel)
            {
                const int bo = org(basel);
                const int bd = dest(basel);
                auto validE = [&](int e) { return ccw(dest(e), bd, bo); };

                int lc = onext(sym(basel));
                bool lv = validE(lc);
                if (lv)
                {
                    bool deleted = false;
                    while (inCircle(bd, bo, dest(lc), dest(onext(lc))))
                    {
                        const int t = onext(lc);
                        deleteEdge(lc);
                        lc = t;
                        deleted = true;
                    }
                    if (deleted) lv = validE(lc);
                }
                return { lc, lv };
            }

            std::pair<int, bool> merge_scan_right(int basel)
            {
                const int bo = org(basel);
                const int bd = dest(basel);
                auto validE = [&](int e) { return ccw(dest(e), bd, bo); };

                int rc = oprev(basel);
                bool rv = validE(rc);
                if (rv)
                {
                    bool deleted = false;
                    while (inCircle(bd, bo, dest(rc), dest(oprev(rc))))
                    {
                        const int t = oprev(rc);
                        deleteEdge(rc);
                        rc = t;
                        deleted = true;
                    }
                    if (deleted) rv = validE(rc);
                }
                return { rc, rv };
            }

            std::pair<int, int> merge(int ldo, int ldi, int rdi, int rdo, int depth)
            {
                for (;;)
                {
                    if (leftOf(org(rdi), ldi)) ldi = lnext(ldi);
                    else if (rightOf(org(ldi), rdi)) rdi = rprev(rdi);
                    else break;
                }
                int basel = connect(sym(rdi), ldi);
                if (org(ldi) == org(ldo)) ldo = sym(basel);
                if (org(rdi) == org(rdo)) rdo = basel;

                const bool concurrent = policy_.split_merge(depth);

                for (;;)
                {
                    const auto scans = policy_.fork_join(concurrent,
                        [this, basel]() { return merge_scan_left(basel); },
                        [this, basel]() { return merge_scan_right(basel); });
                    const int  lc = scans.first.first;
                    const bool lv = scans.first.second;
                    const int  rc = scans.second.first;
                    const bool rv = scans.second.second;

                    if (!lv && !rv) break;

                    const int bo = org(basel);
                    const int bd = dest(basel);
                    if (!lv || (rv && inCircle(dest(lc), bd, bo, dest(rc))))
                        basel = connect(rc, sym(basel));
                    else
                        basel = connect(sym(basel), sym(lc));
                }
                return { ldo, rdo };
            }

            int hullWalk(int e, int frame, bool wantMax) const
            {
                auto better = [&](int u, int v)
                    {
                        return wantMax ? lessF(V[v], V[u], frame) : lessF(V[u], V[v], frame);
                    };
                for (;;)
                {
                    const int nx = rprev(e);
                    const int pv = sym(oprev(e));
                    if (better(org(nx), org(e))) e = nx;
                    else if (better(org(pv), org(e))) e = pv;
                    else return e;
                }
            }

            std::pair<int, int> build(int lo, int hi, int axis, int outAxis, bool presorted, int depth)
            {
                const int n = hi - lo;
                if (n <= 3) return leaf(lo, hi, outAxis);

                const int mid = lo + n / 2;
                if (!presorted)
                {
                    if (axis == 0)
                    {
                        std::nth_element(V.begin() + lo, V.begin() + mid, V.begin() + hi,
                            [](const Vtx& a, const Vtx& b) {
                                return a.x < b.x || (a.x == b.x && a.y < b.y);
                            });
                    }
                    else
                    {
                        std::nth_element(V.begin() + lo, V.begin() + mid, V.begin() + hi,
                            [](const Vtx& a, const Vtx& b) {
                                return a.y < b.y || (a.y == b.y && a.x > b.x);
                            });
                    }
                }
                const int childAxis = 1 - axis;

                const auto halves = policy_.fork_join(policy_.split_build(depth, n),
                    [this, lo, mid, childAxis, axis, depth]() {
                        return build(lo, mid, childAxis, axis, false, depth + 1);
                    },
                    [this, mid, hi, childAxis, axis, depth]() {
                        return build(mid, hi, childAxis, axis, false, depth + 1);
                    });
                const std::pair<int, int>& L = halves.first;
                const std::pair<int, int>& R = halves.second;

                auto res = merge(L.first, L.second, R.first, R.second, depth);

                if (outAxis != axis)
                {
                    const int le = hullWalk(res.first, outAxis, false);
                    const int reCcw = hullWalk(onext(res.second), outAxis, true);
                    res.first = le;
                    res.second = oprev(reCcw);
                }
                return res;
            }
        };

        // The part of a batch triangulation that does not depend on how the
        // builder runs. Fills V with the distinct points of pts, sorted by
        // (x, y); each carries the index of its lowest-numbered occurrence in
        // pts. Records in res.rep, for every input index, the lowest index with
        // the same coordinates, and sets res.unique_points.
        inline void load_unique_points(const Point* pts, size_t sz, Triangulation& res, std::vector<Vtx>& V)
        {
            const int n = static_cast<int>(sz);

            res.rep.assign(n, 0);
            V.reserve(n);
            for (int i = 0; i < n; ++i)
            {
                Vtx v;
                v.x = pts[i].x;
                v.y = pts[i].y;
                v.id = static_cast<size_t>(i);
                V.push_back(v);
            }
            std::iota(res.rep.begin(), res.rep.end(), static_cast<size_t>(0));
            std::sort(V.begin(), V.end(), [](const Vtx& a, const Vtx& b)
                {
                    if (a.x != b.x) return a.x < b.x;
                    if (a.y != b.y) return a.y < b.y;
                    return a.id < b.id;
                });
            size_t m = 0;
            for (size_t i = 0; i < V.size(); ++i)
            {
                if (m > 0 && V[m - 1].x == V[i].x && V[m - 1].y == V[i].y)
                {
                    res.rep[V[i].id] = V[m - 1].id;
                    continue;
                }
                if (m != i) V[m] = V[i];
                ++m;
            }
            V.erase(V.begin() + m, V.end());
            res.unique_points = static_cast<int>(V.size());
        }
    }
}

#endif