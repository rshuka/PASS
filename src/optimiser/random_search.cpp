#include "pass_bits/optimiser/random_search.hpp"

pass::optimise_result pass::random_search::optimise(
    const pass::problem &problem)
{
  pass::optimise_result result(problem.dimension(), acceptable_fitness_value);

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
  } // Termintation criteria
  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved());

  return result;
}
