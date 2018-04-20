#include "pass_bits/helper/evaluation_time_stall.hpp"
#include <thread> // sleep_for

pass::evaluation_time_stall::evaluation_time_stall(const pass::problem &wrapped_problem)
    : problem(wrapped_problem.lower_bounds, wrapped_problem.upper_bounds, wrapped_problem.name),
      wrapped_problem(wrapped_problem),
      repetitions(1) {}

double pass::evaluation_time_stall::evaluate(const arma::vec &agent) const
{
  assert(repetitions >= 1 && "`repetititions` must be at least 1");

  double objective_value = 0.0;
  for (arma::uword n = 0; n < repetitions; n++)
  {
    objective_value = wrapped_problem.evaluate(agent);
  }
  return objective_value;
}
