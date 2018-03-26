#include "pass_bits/optimiser/random_search.hpp"

pass::random_search::random_search() noexcept
    : optimiser("Random_Search_Algorithm") {}

pass::optimise_result pass::random_search::optimise(
    const pass::problem &problem)
{
  pass::optimise_result result(problem, acceptable_fitness_value);

  pass::stopwatch stopwatch;
  stopwatch.start();

  do
  {
    ++result.iterations;

    arma::vec agent = problem.normalised_random_agents(1);
    const double fitness_value = problem.evaluate_normalised(agent);

    if (fitness_value <= result.fitness_value)
    {
      result.normalised_agent = agent;
      result.fitness_value = fitness_value;
    }
    result.duration = stopwatch.get_elapsed();

  } // Termintation criteria
  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved());

  return result;
}
