#include "pass_bits/helper/stopwatch.hpp"

pass::stopwatch::stopwatch() noexcept
    : start_time(std::chrono::steady_clock::now()) {}

std::chrono::microseconds pass::stopwatch::get_elapsed() const noexcept
{
    return std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::steady_clock::now() - start_time);
}
