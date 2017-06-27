#pragma once

#include "../problem.hpp"

namespace pass {

/**
 * The Rastrigin function is a common *toy* problem with a very small
 * computational cost, used for testing and benchmarking algorithms. It is named
 * after L. A. Rastrigin and was first published 1974 in *Systems of Extremal
 * Control* as 2-dimensional problem and 1991 generalised for n dimensions by H.
 * Mühlenbein, D. Schomisch and J. Born in *The Parallel Genetic Algorithm as
 * Function Optimiser. Parallel Computing*. Its optimal parameter = (0, ..., 0)
 * and optimal function value = 0.
 *
 *              N ⎛                            ⎞
 *     10 * N + ∑ ⎜p(i)^2 - 10 * cos(2π * p(i))⎟
 *             i=1⎝                            ⎠
 */
class rastrigin_function : public problem {
 public:
  /**
   * Initialises a sphere function with `dimension` dimensions, lower bounds of
   * -5.12 and upper bounds of 5.12.
   */
  rastrigin_function(const arma::uword dimension);

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
