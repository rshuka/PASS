#include "problem.hpp"

// NAN, INFINITY
#include <cmath>

// assert
#include <cassert>

arma::vec pass::problem::bounds_range() const noexcept {
  return upper_bounds - lower_bounds;
}

std::size_t pass::problem::dimension() const noexcept {
  return lower_bounds.n_elem;
}

pass::problem::problem(const std::size_t dimension, const double lower_bound,
                       const double upper_bound)
    : lower_bounds(arma::vec(dimension).fill(lower_bound)),
      upper_bounds(arma::vec(dimension).fill(upper_bound)) {
  assert(lower_bound < upper_bound &&
         "`problem.lower_bounds` must be less than `problem.upper_bounds`");
}

pass::problem::problem(const arma::vec& lower_bounds,
                       const arma::vec& upper_bounds)
    : lower_bounds(lower_bounds), upper_bounds(upper_bounds) {
  assert(upper_bounds.n_elem == dimension() &&
         "`problems.lower_bounds` and `problem.upper_bounds` must have the "
         "same dimension");
  assert(arma::all(lower_bounds < upper_bounds) &&
         "`problem.lower_bounds` must be less than `problem.upper_bounds` in "
         "each dimension");
}
