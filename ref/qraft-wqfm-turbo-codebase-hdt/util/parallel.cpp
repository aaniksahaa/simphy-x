#include "internal.hpp"

// glibc arena controls do not affect jemalloc and would pull glibc's malloc
// object into an otherwise valid full-static jemalloc link.
#if defined(__GLIBC__) && !defined(WQFM_LINKED_JEMALLOC)
#include <malloc.h>
#endif

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <thread>

namespace wqfm::internal {

namespace {
std::unique_ptr<tf::Executor> g_executor;

#if defined(__GLIBC__) && !defined(WQFM_LINKED_JEMALLOC)
int malloc_arena_cap(int workers) {
  const char* env = std::getenv("WQFM_MALLOC_ARENA_MAX");
  if (env != nullptr && *env != '\0') {
    char* end = nullptr;
    const long parsed = std::strtol(env, &end, 10);
    if (end != env) {
      return std::max(1, static_cast<int>(parsed));
    }
  }
  return std::max(1, std::min(workers, 2));
}
#endif
}

/// Initializes the shared executor from the requested thread count.
int configure_threads(const Options& options) {
  int nw = options.threads > 0 ? options.threads
                               : static_cast<int>(std::thread::hardware_concurrency());
  if (nw < 1) nw = 1;
#if defined(__GLIBC__) && !defined(WQFM_LINKED_JEMALLOC)
  mallopt(M_ARENA_MAX, malloc_arena_cap(nw));
  mallopt(M_TRIM_THRESHOLD, 4 * 1024 * 1024);
#endif
  g_executor = std::make_unique<tf::Executor>(nw);
  return nw;
}

/// Returns the shared executor, creating a default one if needed.
tf::Executor& get_executor() {
  if (!g_executor) {
    g_executor = std::make_unique<tf::Executor>();
  }
  return *g_executor;
}

/// Returns the configured executor worker count.
int get_num_workers() {
  return static_cast<int>(get_executor().num_workers());
}

/// Executes an inner taskflow from worker or external threads.
void run_inner(tf::Taskflow& taskflow) {
  auto& exec = get_executor();
  if (exec.this_worker_id() >= 0) {
    exec.corun(taskflow);
  } else {
    exec.run(taskflow).wait();
  }
}

}  // namespace wqfm::internal
