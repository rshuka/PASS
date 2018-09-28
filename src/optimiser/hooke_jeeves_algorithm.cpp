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

  pass::optimise_result result(problem, acceptable_fitness_value);

  result.normalised_agent = problem.normalised_random_agents(1);
  result.fitness_value = problem.evaluate_normalised(result.normalised_agent);
  ++result.iterations;
  ++result.evaluations;

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
      verbose(result.iterations, 2) = result.normalised_agent[0];
    }
    else
    {
      throw std::runtime_error(
          "Please set - maximal_iterations - to a valid number to analyse the behaviour of the algorithm.");
    }
  }

  double stepsize = 1.0;

  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved())
  {
    arma::vec best_neighbour = result.normalised_agent;
    double neighbour_fitness_value = result.fitness_value;

    for (arma::uword n = 0; n < problem.dimension(); ++n)
    {
      const std::array<double, 2> directions{{stepsize, -stepsize}};

      for (double step : directions)
      {
        arma::vec agent = result.normalised_agent;
        agent(n) += step;

        // place into the boundaries
        if (agent(n) < 0)
        {
          agent(n) = 0;
        }
        else if (agent(n) > 1)
        {
          agent(n) = 1;
        }

        double fitness_value = problem.evaluate_normalised(agent);
        ++result.evaluations;

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
      result.normalised_agent = best_neighbour;
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
      verbose(result.iterations, 2) = result.normalised_agent[0];
    }

    result.duration = stopwatch.get_elapsed();
  }
  // Save the file
  if (pass::is_verbose)
  {
    verbose.shed_row(0);
    verbose.save("Verbose_Optimiser_" + name + "_Problem_" + problem.name + "_Dim_" +
                     std::to_string(problem.dimension()) +
                     "_Run_" + std::to_string(pass::global_number_of_runs),
                 arma::raw_ascii);
  }

  return result;
}
