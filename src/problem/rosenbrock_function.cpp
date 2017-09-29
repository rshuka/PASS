#include "../../include/pass_bits/problem/rosenbrock_function.hpp"

// std::accumulate
#include <algorithm>

// std::sqrt, std::cos
#include <cmath>

// assert
#include <cassert>

pass::rosenbrock_function::rosenbrock_function(const arma::uword dimension)
    : problem(dimension, -2.048, 2.048) {}

double pass::rosenbrock_function::evaluate(const arma::vec& parameter) const {
  assert(parameter.n_elem == dimension() &&
         "`parameter` has incompatible dimension");
  return std::inner_product(
      parameter.cbegin(), std::prev(parameter.cend(), 1),
      std::next(parameter.cbegin(), 1), 0.0, std::plus<double>(),
      [](const double element, const double other_element) {
        return 100.0 * std::pow(other_element - std::pow(element, 2.0), 2.0) +
               std::pow(element - 1.0, 2.0);
      });
}
