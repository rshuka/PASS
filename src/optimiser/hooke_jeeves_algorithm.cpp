#include "pass_bits/optimiser/hooke_jeeves_algorithm.hpp"
#include <array>

pass::optimise_result pass::hooke_jeeves_algorithm::optimise(
    const pass::problem &problem)
{
  pass::stopwatch stopwatch;

  pass::optimise_result result(problem.dimension(), acceptable_objective_value);

  result.parameter = problem.random_parameters(1);
  result.objective_value = problem.evaluate(result.parameter);

  // get the initial stepsize from the dimension with the biggest boundary
  double stepsize = arma::max(problem.bounds_range()) / 2;

  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved())
  {
    arma::vec best_neighbour = result.parameter;
    double neighbour_objective_value = result.objective_value;

    for (arma::uword n = 0; n < problem.dimension(); ++n)
    {
      const std::array<double, 2> directions{{stepsize, -stepsize}};

      for (double step : directions)
      {
        arma::vec parameter = result.parameter;
        parameter(n) += step;

        // place into the boundaries
        if (parameter(n) < problem.lower_bounds(n))
        {
          parameter(n) = problem.lower_bounds(n);
        }
        else if (parameter(n) > problem.upper_bounds(n))
        {
          parameter(n) = problem.upper_bounds(n);
        }

        double objective_value = problem.evaluate(parameter);

        if (objective_value < neighbour_objective_value)
        {
          best_neighbour = parameter;
          neighbour_objective_value = objective_value;
        }
      }
    }

    if (result.objective_value >= neighbour_objective_value)
    {
      // if not improving, halve the stepsize
      stepsize /= 2;
    }
    else
    {
      result.parameter = best_neighbour;
      result.objective_value = neighbour_objective_value;
    }

    ++result.iterations;
    result.duration = stopwatch.get_elapsed();
  }

  return result;
}
