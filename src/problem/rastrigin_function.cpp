#include "../../include/pass_bits/problem/rastrigin_function.hpp"

// std::accumulate
#include <algorithm>

// std::acos, std::cos, std::pow
#include <cmath>

// assert
#include <cassert>

pass::rastrigin_function::rastrigin_function(const arma::uword dimension)
    : problem(dimension, -5.12, 5.12) {}

double pass::rastrigin_function::evaluate(const arma::vec& parameter) const {
  assert(parameter.n_elem == dimension() &&
         "`parameter` has incompatible dimension");
  return 10.0 * dimension() +
         std::accumulate(parameter.cbegin(), parameter.cend(), 0.0,
                         [](const double sum, const double element) {
                           return sum + std::pow(element, 2.0) -
                                  10.0 *
                                      std::cos(2.0 * std::acos(-1.0) * element);
                         });
}
