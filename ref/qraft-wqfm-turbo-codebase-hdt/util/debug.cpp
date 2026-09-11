#include "internal.hpp"

#include <atomic>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>

namespace wqfm::internal {
namespace {

/// Shared debug timing and counter state.
struct DebugState {
  std::mutex mutex;
  std::map<std::string, long long> counters;
  std::map<std::string, long long> maxima;
  std::map<std::string, std::chrono::nanoseconds> timers;
  std::map<std::string, long long> timer_calls;
};

/// Returns the process-global debug state.
DebugState& state() {
  static DebugState s;
  return s;
}

/// Parses environment-style boolean values.
bool env_truthy(const char* value) {
  if (value == nullptr || value[0] == '\0') {
    return false;
  }
  const std::string text(value);
  return text != "0" && text != "false" && text != "FALSE" && text != "off" &&
         text != "OFF";
}

/// Formats a duration in seconds for debug output.
std::string seconds_string(std::chrono::nanoseconds elapsed) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(6)
      << std::chrono::duration<double>(elapsed).count();
  return out.str();
}

}  // namespace

/// Global switch for debug timing and counters.
std::atomic<bool> g_debug_timing_enabled(
  env_truthy(std::getenv("WQFM_DEBUG_TIMING")) || env_truthy(std::getenv("WQFM_DEBUG")));

/// Emits one debug line to stderr.
void debug_log_line_impl(std::string_view message) {
  DebugState& s = state();
  std::lock_guard<std::mutex> lock(s.mutex);
  std::cerr << "[wqfm-debug] " << message << '\n';
}

/// Adds to a named debug counter.
void debug_counter_add_impl(std::string_view name, long long delta) {
  DebugState& s = state();
  std::lock_guard<std::mutex> lock(s.mutex);
  s.counters[std::string(name)] += delta;
}

/// Updates a named debug maximum.
void debug_counter_max_impl(std::string_view name, long long value) {
  DebugState& s = state();
  std::lock_guard<std::mutex> lock(s.mutex);
  auto& current = s.maxima[std::string(name)];
  current = std::max(current, value);
}

/// Adds elapsed time to a named debug timer.
void debug_timer_add_impl(std::string_view name, std::chrono::nanoseconds elapsed) {
  DebugState& s = state();
  std::lock_guard<std::mutex> lock(s.mutex);
  const std::string key(name);
  s.timers[key] += elapsed;
  ++s.timer_calls[key];
}

/// Emits all debug timers, counters, and maxima.
void debug_report_impl() {
  DebugState snapshot;
  {
    DebugState& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    snapshot.counters = s.counters;
    snapshot.maxima = s.maxima;
    snapshot.timers = s.timers;
    snapshot.timer_calls = s.timer_calls;
  }

  std::cerr << "[wqfm-debug] ==== timing summary ====\n";
  for (const auto& [name, elapsed] : snapshot.timers) {
    const long long calls = snapshot.timer_calls[name];
    std::chrono::nanoseconds avg(0);
    if (calls > 0) {
      avg = std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed / calls);
    }
    std::cerr << "[wqfm-debug] timer " << name << " calls=" << calls
              << " total_s=" << seconds_string(elapsed)
              << " avg_s=" << seconds_string(avg) << '\n';
  }
  for (const auto& [name, value] : snapshot.counters) {
    std::cerr << "[wqfm-debug] counter " << name << " value=" << value << '\n';
  }
  for (const auto& [name, value] : snapshot.maxima) {
    std::cerr << "[wqfm-debug] max " << name << " value=" << value << '\n';
  }
}

}  // namespace wqfm::internal
