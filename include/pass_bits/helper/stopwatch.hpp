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
  /**
   * Creat time of the object and start it.
   */
  void start() noexcept;

  /**
   * The time in nanoseconds (10^-9) since this object was created.
   */
  std::chrono::nanoseconds get_elapsed() const noexcept;

private:
  std::chrono::steady_clock::time_point start_time;
};
} // namespace pass
