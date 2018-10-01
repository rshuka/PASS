#include "pass_bits/optimiser.hpp"
#include <cmath> // NAN, INFINITY

pass::optimise_result::optimise_result(const pass::problem &problem, const double acceptable_fitness_value) noexcept
    : normalised_agent(arma::vec(problem.dimension()).fill(std::numeric_limits<double>::quiet_NaN())),
      fitness_value(std::numeric_limits<double>::infinity()),
      acceptable_fitness_value(acceptable_fitness_value),
      problem(problem),
      iterations(0),
      evaluations(0),
      duration(std::chrono::nanoseconds(0)) {}

bool pass::optimise_result::solved() const
{
  return fitness_value <= acceptable_fitness_value;
}

arma::vec pass::optimise_result::agent() const
{
  return normalised_agent % problem.bounds_range() + problem.lower_bounds;
}

pass::optimiser::optimiser(const std::string name)
    : acceptable_fitness_value(-std::numeric_limits<double>::infinity()),
      maximal_iterations(std::numeric_limits<arma::uword>::max()),
      maximal_evaluations(std::numeric_limits<arma::uword>::max()),
      maximal_duration(std::chrono::system_clock::duration::max().count()),
      name(name)
{
  assert(maximal_iterations > 0 &&
         "`maximal_iterations` should be greater than 0");
  assert(maximal_evaluations > 0 &&
         "`maximal_evaluations` should be greater than 0");
  assert(maximal_duration.count() > 0.0 &&
         "`maximal_duration` should be greater than 0");
  assert(name.empty() == false &&
         "`name` should should not be empty");
}
