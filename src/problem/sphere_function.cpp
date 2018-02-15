#include "pass_bits/problem/sphere_function.hpp"

pass::sphere_function::sphere_function(const arma::uword dimension)
    : problem(dimension, -5.12, 5.12) {}

double pass::sphere_function::evaluate(const arma::vec &parameter) const
{
  assert(parameter.n_elem == dimension() &&
         "`parameter` has incompatible dimension");
  return std::accumulate(parameter.begin(), parameter.end(), 0.0,
                         [](const double sum, const double element) {
                           return sum + std::pow(element, 2.0);
                         });
}