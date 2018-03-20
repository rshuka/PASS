#include "pass_bits/problem/sum_of_different_powers_function.hpp"

pass::sum_of_different_powers_function::sum_of_different_powers_function(
    const arma::uword dimension)
    : problem(dimension, -1.0, 1.0, "Sum_Of_Different_Powers_Function") {}

double pass::sum_of_different_powers_function::evaluate(
    const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");

  double sum = 0.0;
  for (std::size_t n = 0; n < dimension(); ++n)
  {
    sum += std::pow(std::fabs(agent(n)), n + 2);
  }
  return sum;
}
