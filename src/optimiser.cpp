#include "../include/pass_bits/optimiser.hpp"

// NAN, INFINITY
#include <cmath>

pass::optimise_result::optimise_result(const arma::uword dimension) noexcept
    : parameter(
          arma::vec(dimension).fill(std::numeric_limits<double>::quiet_NaN())),
      objective_value(std::numeric_limits<double>::infinity()),
      evaluations(0),
      iterations(0),
      duration(std::chrono::nanoseconds(0)) {}

pass::optimiser::optimiser() noexcept
    : acceptable_objective_value(-std::numeric_limits<double>::infinity()),
      maximal_evaluations(std::numeric_limits<arma::uword>::max()),
      maximal_iterations(std::numeric_limits<arma::uword>::max()),
      maximal_duration(std::chrono::minutes(1)) {}
