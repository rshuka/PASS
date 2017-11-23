#pragma once

#include <chrono>

namespace pass {
/**
 * A stopwatch measures the time since its creation.
 */
class stopwatch {
 public:
  stopwatch() noexcept;

  /**
   * The time in nanoseconds since this object was created.
   */
  std::chrono::nanoseconds get_elapsed() const noexcept;

 private:
  const std::chrono::steady_clock::time_point start_time;
};
}  // namespace pass
