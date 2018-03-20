#include "pass_bits/optimiser.hpp"
#include <cmath> // NAN, INFINITY

pass::optimise_result::optimise_result(const arma::uword dimension, const double acceptable_fitness_value) noexcept
    : agent(arma::vec(dimension).fill(std::numeric_limits<double>::quiet_NaN())),
      fitness_value(std::numeric_limits<double>::infinity()),
      acceptable_fitness_value(acceptable_fitness_value),
      iterations(0),
      duration(std::chrono::nanoseconds(0)) {}

bool pass::optimise_result::solved() const
{
  return fitness_value <= acceptable_fitness_value;
}

pass::optimiser::optimiser(const std::string name)
    : acceptable_fitness_value(-std::numeric_limits<double>::infinity()),
      maximal_iterations(std::numeric_limits<arma::uword>::max()),
      maximal_duration(std::chrono::system_clock::duration::max()),
      name(name)
{
  assert(maximal_iterations > 0 &&
         "`maximal_iterations` should be greater than 0");
  assert(maximal_duration > 0 &&
         "`maximal_iterations` should be greater than 0");
  assert(name != name.empty() &&
         "`name` should should not be empty");
}
