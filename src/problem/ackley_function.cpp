#include "pass_bits/problem/ackley_function.hpp"

pass::ackley_function::ackley_function(const arma::uword dimension)
    : problem(dimension, -32.768, 32.768) {}

double pass::ackley_function::evaluate(const arma::vec &parameter) const
{
  assert(parameter.n_elem == dimension() &&
         "`parameter` has incompatible dimension");
  return -20.0 * std::exp(-0.2 * std::sqrt(1.0 / dimension() *
                                           arma::sum(parameter % parameter))) -
         std::exp(1.0 / dimension() *
                  std::accumulate(
                      parameter.cbegin(), parameter.cend(), 0.0,
                      [](const double sum, const double element) {
                        return sum + std::cos(2.0 * arma::datum::pi * element);
                      })) +
         20.0 + std::exp(1.0);
}
