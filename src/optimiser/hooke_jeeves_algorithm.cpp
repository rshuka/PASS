#include "pass_bits/optimiser/hooke_jeeves_algorithm.hpp"

#include <array>
#include <cassert>

pass::hooke_jeeves_algorithm::hooke_jeeves_algorithm() noexcept
    : initial_stepsize(1), stepsize_decrease(2) {}

pass::optimise_result pass::hooke_jeeves_algorithm::optimise(
    const pass::problem& problem) {
  assert(stepsize_decrease > 1.0 && "");

  pass::optimise_result result(problem.dimension());
  auto start_time = std::chrono::steady_clock::now();

  result.parameter = problem.random_parameters(1);
  result.objective_value = problem.evaluate(result.parameter);
  ++result.evaluations;

  double stepsize = initial_stepsize;

  while (result.duration <= maximal_duration &&
         result.evaluations <= maximal_evaluations &&
         result.objective_value > acceptable_objective_value) {
    bool is_improving = false;

    for (std::size_t n = 0; n < problem.dimension(); ++n) {
      const std::array<double, 2> directions{{stepsize, -2 * stepsize}};
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

        if (objective_value < result.objective_value) {
          is_improving = true;

          result.parameter = parameter;
          result.objective_value = objective_value;

          if (result.objective_value <= this->acceptable_objective_value) {
            result.solved = true;
            return result;
          }
        }

        if (result.evaluations >= this->maximal_evaluations) {
          result.solved = result.objective_value <= acceptable_objective_value;
          return result;
        } else if (result.duration >= this->maximal_duration) {
          result.solved = result.objective_value <= acceptable_objective_value;
          return result;
        }
      }
    }
    ++result.iterations;

    if (!is_improving) {
      stepsize /= stepsize_decrease;
    }
  }

  result.solved = result.objective_value <= acceptable_objective_value;
  return result;
}
