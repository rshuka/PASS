#include "pass_bits/optimiser/hooke_jeeves_algorithm.hpp"
#include "pass_bits/helper/stopwatch.hpp"

#include <array>
#include <cassert>

pass::hooke_jeeves_algorithm::hooke_jeeves_algorithm() noexcept
    : initial_stepsize(0), stepsize_decrease(2) {}

pass::optimise_result pass::hooke_jeeves_algorithm::optimise(
    const pass::problem& problem) {
  assert(initial_stepsize > 0.0 && "initial_stepsize must not be 0");
  assert(stepsize_decrease > 1.0 && "stepsize_decrease must be greater than 1");

  pass::stopwatch stopwatch;

  pass::optimise_result result(problem.dimension(), acceptable_objective_value,
                               2 * problem.dimension());

  result.parameter = problem.random_parameters(1);
  result.objective_value = problem.evaluate(result.parameter);

  double stepsize = initial_stepsize;

  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved()) {
    arma::vec best_neighbour = result.parameter;
    double neighbour_objective_value = result.objective_value;

    for (std::size_t n = 0; n < problem.dimension(); ++n) {
      const std::array<double, 2> directions{{stepsize, -stepsize}};
      for (double step : directions) {
        arma::vec parameter = result.parameter;
        parameter(n) += step;
        if (parameter(n) < problem.lower_bounds(n)) {
          parameter(n) = problem.lower_bounds(n);
        } else if (parameter(n) > problem.upper_bounds(n)) {
          parameter(n) = problem.upper_bounds(n);
        }

        double objective_value = problem.evaluate(parameter);
        if (objective_value < neighbour_objective_value) {
          best_neighbour = parameter;
          neighbour_objective_value = objective_value;
        }
      }
    }

    if (result.objective_value == neighbour_objective_value) {
      // not improving
      stepsize /= stepsize_decrease;
    } else {
      result.parameter = best_neighbour;
      result.objective_value = neighbour_objective_value;
    }

    ++result.iterations;
    result.duration stopwatch.get_elapsed();
  }

  return result;
}
