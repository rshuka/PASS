#include "pass_bits/optimiser/random_search.hpp"
#include "pass_bits/helper/stopwatch.hpp"

pass::optimise_result pass::random_search::optimise(
    const pass::problem &problem)
{
  pass::optimise_result result(problem.dimension(), acceptable_objective_value);
  pass::stopwatch stopwatch;

  do
  {
    ++result.iterations;

    arma::vec parameter = problem.random_parameters(1);
    const double objective_value = problem.evaluate(parameter);

    if (objective_value <= result.objective_value)
    {
      result.parameter = parameter;
      result.objective_value = objective_value;
    }

  } // Termintation criteria
  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved());

  result.duration = stopwatch.get_elapsed();

  return result;
}
