#pragma once

#include "../problem.hpp"

namespace pass {

/**
 * The Rosenbrock function is a common *toy* problem with a very small
 * computational cost, used for testing and benchmarking algorithms. It is named
 * after Howard H. Rosenbrock and was first published 1960 in *An automatic
 * method for finding the greatest or least value of a function. The Computer
 * Journal*. Its optimal parameter = (1, ..., 1) and optimal function value = 0.
 *
 *         ⎛      ⎛                     ⎞     ⎛          ⎞   ⎞
 *     N-1 ⎜      ⎜                     ⎟^2   ⎜          ⎟^2 ⎟
 *      ∑  ⎜100 * ⎜ p(i + 1) - (p(i))^2 ⎟   + ⎜ p(i) - 1 ⎟   ⎟
 *     i=1 ⎜      ⎜                     ⎟     ⎜          ⎟   ⎟
 *         ⎝      ⎝                     ⎠     ⎝          ⎠   ⎠
 */
class rosenbrock_function : public problem {
 public:
  /**
   * Initialises a sphere function with `dimension` dimensions, lower bounds of
   * -2.048 and upper bounds of 2.048.
   */
  rosenbrock_function(const arma::uword dimension);

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
