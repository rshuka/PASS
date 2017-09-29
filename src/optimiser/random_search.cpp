#include "../../include/pass_bits/optimiser/random_search.hpp"

#include <cassert>

pass::optimise_result pass::random_search::optimise(
    const pass::problem& problem) {
  pass::optimise_result result(problem.dimension());
  auto start_time = std::chrono::steady_clock::now();

  do {
    arma::vec parameter = problem.random_parameters(1).col(0);
    const double objective_value = problem.evaluate(parameter);

    if (objective_value <= result.objective_value) {
      result.parameter = parameter;
      result.objective_value = objective_value;
    }

    result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start_time);
    ++result.evaluations;
    ++result.iterations;
  }  // Termintation criteria
  while (result.duration <= maximal_duration &&
         result.evaluations <= maximal_evaluations &&
         result.objective_value > acceptable_objective_value);

  return result;
}
