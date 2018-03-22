#include "pass_bits/problem/optimisation_benchmark/styblinski_tang_function.hpp"

pass::styblinski_tang_function::styblinski_tang_function(const arma::uword dimension)
    : problem(dimension, -5.0, 5.0, "Styblinski_Tang_Function") {}

double pass::styblinski_tang_function::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");
  return 0.5 *
         std::accumulate(agent.cbegin(), agent.cend(), 0.0,
                         [](const double sum, const double element) {
                           return sum + std::pow(element, 4) - 16 * std::pow(element, 2) + 5 * element;
                         });
}
