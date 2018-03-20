#include "pass_bits/optimiser/random_search.hpp"

pass::random_search::random_search() noexcept
    : optimiser("Random_Search_Algorithm") {}

pass::optimise_result pass::random_search::optimise(
    const pass::problem &problem)
{
  pass::optimise_result result(problem.dimension(), acceptable_fitness_value);

  // Variables used to analyse the behavior of a particle
  arma::mat verbose(maximal_iterations + 1, 3);

  pass::stopwatch stopwatch;
  stopwatch.start();

  do
  {
    ++result.iterations;

    arma::vec agent = problem.random_agents(1);
    const double fitness_value = problem.evaluate(agent);

    if (fitness_value <= result.fitness_value)
    {
      result.agent = agent;
      result.fitness_value = fitness_value;
    }
    result.duration = stopwatch.get_elapsed();

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

  } // Termintation criteria
  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved());

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
