#include "optimiser.hpp"

// NAN, INFINITY
#include <cmath>

// assert
#include <cassert>

pass::optimise_result::optimise_result(const std::size_t dimension) noexcept
    : parameter(
          arma::vec(dimension).fill(std::numeric_limits<double>::quiet_NaN())),
      objective_value(std::numeric_limits<double>::infinity()),
      evaluations(0),
      duration(std::chrono::nanoseconds(0)) {}

pass::optimiser::optimiser() noexcept
    : acceptable_objective_value(-std::numeric_limits<double>::infinity()),
      maximal_evaluations(std::numeric_limits<std::size_t>::max()),
      maximal_duration(std::chrono::seconds(1)) {}
