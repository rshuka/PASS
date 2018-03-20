#include "pass_bits/optimiser/hooke_jeeves_algorithm.hpp"
#include <array>

pass::hooke_jeeves_algorithm::hooke_jeeves_algorithm() noexcept
    : optimiser("Hooke_Jeeves_Algorithm") {}

pass::optimise_result pass::hooke_jeeves_algorithm::optimise(
    const pass::problem &problem)
{
  // Variables used to analyse the behavior of a particle
  arma::mat verbose(maximal_iterations + 1, 3);

  pass::stopwatch stopwatch;
  stopwatch.start();

  pass::optimise_result result(problem.dimension(), acceptable_fitness_value);

  result.agent = problem.random_agents(1);
  result.fitness_value = problem.evaluate(result.agent);
  ++result.iterations;

  /**
    * +------------+---------------+----------+
    * | Iterations | Fitness Value | Position |
    * +------------+---------------+----------+
    * Each Dimension is independent. So, the analysis can be performed
    * on just one dimension
    */
  if (pass::is_verbose)
  {
    if (maximal_iterations != std::numeric_limits<arma::uword>::max() && maximal_iterations > 0)
    {
      verbose(result.iterations, 0) = result.iterations;
      verbose(result.iterations, 1) = result.fitness_value;
      verbose(result.iterations, 2) = result.agent[0];
    }
    else
    {
      throw std::runtime_error(
          "Please set - maximal_iterations - to a valid number to analyse the behaviour of the algorithm.");
    }
  }

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

    /**
     * +------------+--------------+----------+
     * | Iterations |Fitness Value | Position |
     * +------------+--------------+----------+
     * Each Dimension is independent. So, the analysis can be performed
     * on just one dimension
     */
    if (pass::is_verbose)
    {
      verbose(result.iterations, 0) = result.iterations;
      verbose(result.iterations, 1) = result.fitness_value;
      verbose(result.iterations, 2) = result.agent[0];
    }

    result.duration = stopwatch.get_elapsed();
  }
  // Save the file
  if (pass::is_verbose)
  {
    verbose.shed_row(0);
    verbose.save("Verbose_" + name + "_Problem_" + problem.name + "_Dim_" +
                     std::to_string(problem.dimension()) +
                     "_Run_" + std::to_string(pass::number_of_runs),
                 arma::raw_ascii);
  }

  return result;
}
