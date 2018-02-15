#include "pass_bits/problem/sum_of_different_powers.hpp"

pass::sum_of_different_powers::sum_of_different_powers(
    const arma::uword dimension)
    : problem(dimension, -1.0, 1.0) {}

double pass::sum_of_different_powers::evaluate(
    const arma::vec &parameter) const
{
  assert(parameter.n_elem == dimension() &&
         "`parameter` has incompatible dimension");

  double sum = 0.0;
  for (std::size_t n = 0; n < dimension(); ++n)
  {
    sum += std::pow(std::fabs(parameter(n)), n + 2);
  }
  return sum;
}