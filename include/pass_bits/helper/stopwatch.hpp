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
   * The time in microseconds (10^-6) since this object was created.
   */
  std::chrono::microseconds get_elapsed() const noexcept;

private:
  std::chrono::steady_clock::time_point start_time;
};
} // namespace pass
