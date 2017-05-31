#include "../../include/pass_bits/problem/ackley_function.hpp"

// std::accumulate
#include <algorithm>

// std::sqrt, std::cos
#include <cmath>

// assert
#include <cassert>

pass::ackley_function::ackley_function(const arma::uword dimension)
    : problem(dimension, -32.768, 32.768) {}

double pass::ackley_function::evaluate(const arma::vec& parameter) const {
  assert(parameter.n_elem == dimension() &&
         "`parameter` has incompatible dimension");
  return 20.0 *
             (-std::exp(-0.2 * arma::norm(parameter) / std::sqrt(dimension())) +
              1.0) -
         std::exp(arma::sum(std::cos(2 * arma::datum::pi) * parameter) /
                  dimension()) +
         std::exp(1.0);
}
