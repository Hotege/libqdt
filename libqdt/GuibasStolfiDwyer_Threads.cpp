#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <future>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "qdt.h"
#include "builder.hxx"
#include "robust.hxx"
#include "wrapper.hxx"

namespace parallel
{
    constexpr unsigned MAX_WORKERS = 64;

    thread_local int tls_slot = 0;

    // Thrown when the parallel builder cannot continue (its pre-allocated edge
    // storage is exhausted, or a thread slot is out of range). It never leaves
    // this file: delaunay_triangulate() catches it and falls back to the serial
    // implementation.
    struct BuilderOverflow : std::runtime_error
    {
        explicit BuilderOverflow(const char* what) : std::runtime_error(what) {}
    };

    class ThreadPool
    {
    public:
        explicit ThreadPool(size_t n) : stop_(false)
        {
            if (n == 0) n = 1;
            try
            {
                for (size_t i = 0; i < n; ++i)
                    workers_.emplace_back([this, i] {
                    tls_slot = static_cast<int>(i) + 1;
                    worker_loop();
                        });
            }
            catch (...)
            {
                // The destructor does not run for a half-constructed object,
                // and destroying a joinable std::thread terminates the process.
                // Stop and join the workers that did start, then rethrow.
                shutdown_();
                throw;
            }
        }

        size_t size() const { return workers_.size(); }

        ~ThreadPool() { shutdown_(); }

        template <class F>
        auto enqueue(F&& f) -> std::future<decltype(f())>
        {
            using R = decltype(f());
            auto task = std::make_shared<std::packaged_task<R()>>(
                std::forward<F>(f));
            std::future<R> fut = task->get_future();
            {
                std::unique_lock<std::mutex> lock(mtx_);
                tasks_.emplace([task]() { (*task)(); });
            }
            cv_.notify_one();
            return fut;
        }

        bool try_run_one()
        {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(mtx_);
                if (tasks_.empty()) return false;
                task = std::move(tasks_.front());
                tasks_.pop();
            }
            task();
            return true;
        }

        template <class T>
        T wait_for(std::future<T>& fut)
        {
            for (;;)
            {
                if (fut.wait_for(std::chrono::seconds(0)) ==
                    std::future_status::ready)
                    return fut.get();
                if (!try_run_one())
                    return fut.get();
            }
        }

    private:
        void shutdown_()
        {
            {
                std::unique_lock<std::mutex> lock(mtx_);
                stop_ = true;
            }
            cv_.notify_all();
            for (auto& t : workers_)
            {
                if (t.joinable()) t.join();
            }
        }

        void worker_loop()
        {
            for (;;)
            {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(mtx_);
                    cv_.wait(lock, [this] { return stop_ || !tasks_.empty(); });
                    if (stop_ && tasks_.empty()) return;
                    task = std::move(tasks_.front());
                    tasks_.pop();
                }
                task();
            }
        }

        std::vector<std::thread> workers_;
        std::queue<std::function<void()>> tasks_;
        std::mutex mtx_;
        std::condition_variable cv_;
        bool stop_;
    };

    using wrapper::Point;

    // Edge storage and concurrency of the parallel builder. See builder.hxx
    // for what a policy has to provide.
    //
    // The quad-edge array is allocated in full up front (5 records per point
    // plus 64 Ki), and the threads take it in blocks of CHUNK_SIZE records, so
    // that they rarely touch the shared counter. Every thread has its own
    // block and its own free list (slot 0: the calling thread, slots 1..N: the
    // pool workers). If the array runs out, BuilderOverflow is thrown and the
    // caller falls back to the serial implementation.
    class ParallelPolicy
    {
    public:
        static constexpr int CHUNK_SIZE = 65536;
        // slot 0 = calling thread, slots 1..MAX_WORKERS = pool workers
        static constexpr int MAX_POOLS = static_cast<int>(MAX_WORKERS) + 1;
        static constexpr int MAX_PARALLEL_DEPTH = 3;
        static constexpr int PARALLEL_MIN_SIZE = 524288;
        static constexpr int MERGE_PARALLEL_DEPTH = 1;

        // The pool is attached after construction so that the caller can
        // declare the builder before the pool: the pool (which joins its
        // workers) is then destroyed first, and no task can touch a
        // destroyed builder when an exception unwinds the caller.
        void set_pool(ThreadPool* pool) { pool_ = pool; }

        void init(std::vector<wrapper::detail::QE>& qe, size_t nverts)
        {
            qe.resize(5 * static_cast<int>(nverts) + 64 * 1024);
        }

        int alloc(std::vector<wrapper::detail::QE>& qe)
        {
            Pool& p = my_pool();
            if (p.free_head >= 0)
            {
                const int k = p.free_head;
                p.free_head = qe[k].next[0];
                return k;
            }
            if (p.next_slot >= p.block_end)
            {
                const int s = g_next.fetch_add(CHUNK_SIZE, std::memory_order_acq_rel);
                if (s + CHUNK_SIZE > static_cast<int>(qe.size()))
                    throw BuilderOverflow("QDT: parallel edge storage exhausted");
                for (int i = 0; i < CHUNK_SIZE; ++i)
                {
                    qe[s + i].vtx[0] = -1;
                    qe[s + i].vtx[1] = -1;
                }
                p.next_slot = s;
                p.block_end = s + CHUNK_SIZE;
            }
            return p.next_slot++;
        }

        void release(std::vector<wrapper::detail::QE>& qe, int k)
        {
            Pool& p = my_pool();
            qe[k].vtx[0] = -1;
            qe[k].next[0] = p.free_head;
            p.free_head = k;
        }

        int edge_count(const std::vector<wrapper::detail::QE>&) const
        {
            return g_next.load(std::memory_order_acquire);
        }

        bool split_build(int depth, int n) const
        {
            return pool_ != nullptr && depth < MAX_PARALLEL_DEPTH && n >= PARALLEL_MIN_SIZE;
        }

        bool split_merge(int depth) const
        {
            return pool_ != nullptr && depth < MERGE_PARALLEL_DEPTH;
        }

        template <class A, class B>
        auto fork_join(bool concurrent, A&& first, B&& second)
            -> std::pair<decltype(first()), decltype(second())>
        {
            typedef std::pair<decltype(first()), decltype(second())> Result;
            if (!concurrent)
            {
                auto x = first();
                auto y = second();
                return Result(x, y);
            }
            auto fut = pool_->enqueue(std::forward<A>(first));
            auto y = second();
            auto x = pool_->wait_for(fut);
            return Result(x, y);
        }

    private:
        struct alignas(64) Pool
        {
            int next_slot = 0;
            int block_end = 0;
            int free_head = -1;
            char pad[64 - 3 * sizeof(int)];
        };
        static_assert(sizeof(Pool) == 64, "Pool must fit exactly one cache line");

        std::array<Pool, MAX_POOLS> pools_{};
        std::atomic<int> g_next{ 0 };
        ThreadPool* pool_ = nullptr;

        Pool& my_pool()
        {
            const int slot = tls_slot;
            if (slot < 0 || slot >= MAX_POOLS)
                throw BuilderOverflow("QDT: thread slot out of range");
            return pools_[slot];
        }
    };

    wrapper::Triangulation delaunay_triangulate(const Point* pts, size_t sz)
    {
        using Builder = wrapper::detail::Builder<ParallelPolicy>;

        wrapper::Triangulation res;
        std::vector<wrapper::detail::Vtx> V;
        wrapper::detail::load_unique_points(pts, sz, res, V);

        if (V.size() >= 3)
        {
            unsigned nthreads = std::thread::hardware_concurrency();
            if (nthreads == 0) nthreads = 1;
            if (nthreads > MAX_WORKERS) nthreads = MAX_WORKERS;

            if (nthreads > 1 && V.size() >= 200000)
            {
                bool fall_back = false;
                try
                {
                    // The builder is declared before the pool on purpose: after
                    // an exception, tasks that are still queued or running must
                    // not outlive it (see ParallelPolicy::set_pool).
                    Builder b(V);
                    ThreadPool pool(nthreads);
                    b.policy().set_pool(&pool);
                    b.run();
                    b.extract(res.triangles);
                }
                catch (const BuilderOverflow&) { fall_back = true; }
                catch (const std::system_error&) { fall_back = true; } // could not start threads
                if (fall_back)
                {
                    // V was reordered by the failed attempt; the serial
                    // implementation starts again from the caller's points.
                    return wrapper::delaunay_triangulate(pts, sz);
                }
            }
            else
            {
                Builder b(V);
                b.run();
                b.extract(res.triangles);
            }
        }

        res.input_points.resize(sz);
        for (size_t i = 0; i < sz; ++i)
            res.input_points[i] = pts[i];

        return res;
    }
}