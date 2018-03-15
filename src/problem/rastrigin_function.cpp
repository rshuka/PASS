#include "pass_bits/problem/rastrigin_function.hpp"

pass::rastrigin_function::rastrigin_function(const arma::uword dimension)
    : problem(dimension, -5.12, 5.12) {}

double pass::rastrigin_function::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");
  return 10.0 * dimension() +
         std::accumulate(agent.cbegin(), agent.cend(), 0.0,
                         [](const double sum, const double element) {
                           return sum + std::pow(element, 2.0) -
                                  10.0 *
                                      std::cos(2.0 * arma::datum::pi * element);
                         });
}
