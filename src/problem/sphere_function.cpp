#include "sphere_function.hpp"

// std::accumulate
#include <algorithm>

// std::pow
#include <cmath>

// assert
#include <cassert>

pass::sphere_function::sphere_function(const std::size_t dimension)
    : problem(dimension, -5.12, 5.12) {}

double pass::sphere_function::evaluate(const arma::vec& parameter) const {
  assert(parameter.n_elem == dimension() &&
         "`parameter` has incompatible dimension");
  return std::accumulate(parameter.begin(), parameter.end(), 0.0,
                         [](const auto length, const auto element) {
                           return length + std::pow(element, 2.0);
                         });
}
