#include "pass_bits/problem/optimisation_benchmark/griewank_function.hpp"

pass::griewank_function::griewank_function(const arma::uword dimension)
    : problem(dimension, -600.00, 600.00, "Griewank_Function") {}

double pass::griewank_function::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");

  double product = 1.0;
  double sum = 0.0;

  for (arma::uword i = 0; i < agent.n_elem; ++i)
  {
    sum = sum + agent(i) * agent(i);
    product = product * std::cos(agent(i) / std::sqrt(static_cast<double>(i) + 1.0));
  }

  return sum / 4000.0 - product + 1.0;
}
