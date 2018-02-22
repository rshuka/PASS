#pragma once

#include <chrono> // std::chrono

namespace pass
{
/**
 * A stopwatch measures the time since its creation.
 */
class stopwatch
{
public:
  stopwatch() noexcept;

  /**
   * The time in microseconds since this object was created.
   */
  std::chrono::microseconds get_elapsed() const noexcept;

private:
  const std::chrono::steady_clock::time_point start_time;
};
} // namespace pass
