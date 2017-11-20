#include "pass_bits/optimiser/hooke_jeeves_algorithm.hpp"

#include <array>
#include <cassert>

pass::hooke_jeeves_algorithm::hooke_jeeves_algorithm() noexcept
    : initial_stepsize(1), stepsize_decrease(2) {}

pass::optimise_result pass::hooke_jeeves_algorithm::optimise(
    const pass::problem& problem) {
  assert(initial_stepsize > 0.0 && "initial_stepsize must not be 0");
  assert(stepsize_decrease > 1.0 && "stepsize_decrease must be greater than 1");

  pass::optimise_result result(problem.dimension());
  auto start_time = std::chrono::steady_clock::now();

  result.parameter = problem.random_parameters(1);
  result.objective_value = problem.evaluate(result.parameter);
  ++result.evaluations;

  double stepsize = initial_stepsize;

  while (result.duration <= maximal_duration &&
         result.evaluations <= maximal_evaluations &&
         result.objective_value > acceptable_objective_value) {
    arma::vec best_neighbour = result.parameter;
    double neighbour_objective_value = result.objective_value;

    const std::array<double, 2> directions{{stepsize, -stepsize}};

    for (std::size_t n = 0; n < problem.dimension(); ++n) {
      for (double step : directions) {
        arma::vec parameter = result.parameter;
        parameter(n) += step;
        if (parameter(n) < problem.lower_bounds(n)) {
          parameter(n) = problem.lower_bounds(n);
        } else if (parameter(n) > problem.upper_bounds(n)) {
          parameter(n) = problem.upper_bounds(n);
        }

        auto objective_value = problem.evaluate(parameter);
        ++result.evaluations;
        result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - start_time);

        if (objective_value < neighbour_objective_value) {
          if (result.objective_value <= this->acceptable_objective_value) {
            result.parameter = parameter;
            result.objective_value = objective_value;
            result.solved = true;
            return result;
          }

          best_neighbour = parameter;
          neighbour_objective_value = objective_value;
        }

        if (result.evaluations >= this->maximal_evaluations ||
            result.duration >= this->maximal_duration) {
          result.parameter = best_neighbour;
          result.objective_value = neighbour_objective_value;
          return result;
        }
      }
    }
    ++result.iterations;

    if (result.objective_value == neighbour_objective_value) {
      // not improving
      stepsize /= stepsize_decrease;
    } else {
      result.parameter = best_neighbour;
      result.objective_value = neighbour_objective_value;
    }
  }

  result.solved = result.objective_value <= acceptable_objective_value;
  return result;
}
