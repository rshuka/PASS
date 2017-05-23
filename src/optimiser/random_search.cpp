#include "../../include/pass_bits/optimiser/random_search.hpp"

#include <cassert>

pass::optimise_result pass::random_search::optimise(
    const pass::problem& problem, const arma::mat& initial_parameters) {
  assert(initial_parameters.n_elem == 0 ||
         initial_parameters.n_rows == problem.dimension() &&
             "`initial_parameters` must have the same dimension as `problem`, "
             "or be empty");

  auto start_time = std::chrono::steady_clock::now();
  pass::optimise_result result(problem.dimension());
  arma::uword it = 0;

  do {
    arma::vec parameter(problem.dimension());
    if (it < initial_parameters.n_cols) {
      parameter = initial_parameters.col(it++);
    } else {
      parameter.randu();
      parameter %= problem.bounds_range();  // element-wise multiplication
      parameter += problem.lower_bounds;
    }

    const auto objective_value = problem.evaluate(parameter);
    if (objective_value <= result.objective_value) {
      result.parameter = parameter;
      result.objective_value = objective_value;
    }

    ++result.evaluations;
    result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start_time);
  } while (result.duration < maximal_duration &&
           result.evaluations < maximal_evaluations &&
           result.objective_value > acceptable_objective_value);

  return result;
}
