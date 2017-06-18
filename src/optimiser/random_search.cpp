#include "../../include/pass_bits/optimiser/random_search.hpp"

#include <cassert>

pass::optimise_result pass::random_search::optimise(
    const pass::problem& problem) {
  auto start_time = std::chrono::steady_clock::now();
  pass::optimise_result result(problem.dimension());

  do {
    arma::vec parameter(problem.dimension(), arma::fill::randu);
    parameter %= problem.bounds_range();  // element-wise multiplication
    parameter += problem.lower_bounds;

    const auto objective_value = problem.evaluate(problem.random_parameters(1));
    if (objective_value <= result.objective_value) {
      result.parameter = parameter;
      result.objective_value = objective_value;
    }

    ++result.evaluations;
    ++result.iterations;
    result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start_time);
  } while (result.duration < maximal_duration &&
           result.evaluations < maximal_evaluations &&
           result.objective_value > acceptable_objective_value);

  return result;
}
