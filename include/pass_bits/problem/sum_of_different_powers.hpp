#pragma once

// Armadillo
#include <armadillo>

// PASS
#include "pass_bits/problem.hpp"

namespace pass {

/**
 * The sum of different powers function is a common *toy* problem with a very
 * small computational cost, used for testing and benchmarking algorithms. Its
 * optimal parameter = (0, ..., 0) and optimal function value = 0.
 *
 *  N  ⎛                 ⎞
 *  ∑  ⎜||p(i)||^(i + 1) ⎟
 * i=1 ⎝                 ⎠
 */
class sum_of_different_powers : public problem {
 public:
  /**
   * Initialises a different powers function with `dimension` dimensions, lower
   * bounds of -1.0 and upper bounds of 1.0.
   */
  sum_of_different_powers(const arma::uword dimension);

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
