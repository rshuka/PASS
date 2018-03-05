#include "pass_bits/optimiser/hooke_jeeves_algorithm.hpp"
#include <array>

pass::optimise_result pass::hooke_jeeves_algorithm::optimise(
    const pass::problem &problem)
{
  pass::stopwatch stopwatch;
  stopwatch.start();

  pass::optimise_result result(problem.dimension(), acceptable_fitness_value);

  result.agent = problem.random_agents(1);
  result.fitness_value = problem.evaluate(result.agent);

  // get the initial stepsize from the dimension with the biggest boundary
  double stepsize = arma::max(problem.bounds_range()) / 2;

  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved())
  {
    arma::vec best_neighbour = result.agent;
    double neighbour_fitness_value = result.fitness_value;

    for (arma::uword n = 0; n < problem.dimension(); ++n)
    {
      const std::array<double, 2> directions{{stepsize, -stepsize}};

      for (double step : directions)
      {
        arma::vec agent = result.agent;
        agent(n) += step;

        // place into the boundaries
        if (agent(n) < problem.lower_bounds(n))
        {
          agent(n) = problem.lower_bounds(n);
        }
        else if (agent(n) > problem.upper_bounds(n))
        {
          agent(n) = problem.upper_bounds(n);
        }

        double fitness_value = problem.evaluate(agent);

        if (fitness_value < neighbour_fitness_value)
        {
          best_neighbour = agent;
          neighbour_fitness_value = fitness_value;
        }
      }
    }

    if (result.fitness_value >= neighbour_fitness_value)
    {
      // if not improving, halve the stepsize
      stepsize /= 2.0;
    }
    else
    {
      result.agent = best_neighbour;
      result.fitness_value = neighbour_fitness_value;
    }

    ++result.iterations;
    result.duration = stopwatch.get_elapsed();
  }

  return result;
}
