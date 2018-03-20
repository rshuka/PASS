#include "pass_bits/problem/ackley_function.hpp"

pass::ackley_function::ackley_function(const arma::uword dimension)
    : problem(dimension, -32.768, 32.768, "Ackley_Function") {}

double pass::ackley_function::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");
  return -20.0 * std::exp(-0.2 * std::sqrt(1.0 / dimension() *
                                           arma::sum(agent % agent))) -
         std::exp(1.0 / dimension() *
                  std::accumulate(
                      agent.cbegin(), agent.cend(), 0.0,
                      [](const double sum, const double element) {
                        return sum + std::cos(2.0 * arma::datum::pi * element);
                      })) +
         20.0 + std::exp(1.0);
}
