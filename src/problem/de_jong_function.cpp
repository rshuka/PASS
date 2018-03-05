#include "pass_bits/problem/de_jong_function.hpp"

pass::de_jong_function::de_jong_function(const arma::uword dimension)
    : problem(dimension, -5.12, 5.12) {}

double pass::de_jong_function::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");
  return std::accumulate(agent.begin(), agent.end(), 0.0,
                         [](const double sum, const double element) {
                           return sum + std::pow(element, 2.0);
                         });
}