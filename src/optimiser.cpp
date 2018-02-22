#include "pass_bits/optimiser.hpp"
#include <cmath> // NAN, INFINITY

pass::optimise_result::optimise_result(
    const arma::uword dimension,
    const double acceptable_objective_value) noexcept
    : parameter(
          arma::vec(dimension).fill(std::numeric_limits<double>::quiet_NaN())),
      objective_value(std::numeric_limits<double>::infinity()),
      acceptable_objective_value(acceptable_objective_value),
      iterations(0),
      duration(std::chrono::microseconds(0)),
      evaluations_per_iteration(1) {}

pass::optimise_result::optimise_result(
    const arma::uword dimension, const double acceptable_objective_value,
    const arma::uword evaluations_per_iteration) noexcept
    : parameter(
          arma::vec(dimension).fill(std::numeric_limits<double>::quiet_NaN())),
      objective_value(std::numeric_limits<double>::infinity()),
      acceptable_objective_value(acceptable_objective_value),
      iterations(0),
      duration(std::chrono::microseconds(0)),
      evaluations_per_iteration(evaluations_per_iteration) {}

bool pass::optimise_result::solved() const
{
    return objective_value <= acceptable_objective_value;
}

pass::optimiser::optimiser() noexcept
    : acceptable_objective_value(-std::numeric_limits<double>::infinity()),
      maximal_iterations(std::numeric_limits<arma::uword>::max()),
      maximal_duration(std::chrono::minutes(1)) {}