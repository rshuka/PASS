#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * The Rosenbrock function, also referred to as the Valley or Banana function,
 * is a popular test problem for gradient-based optimization algorithms.
 * The function is unimodal, and the global minimum lies in a narrow, parabolic
 * valley. However, even though this valley is easy to find, convergence to the
 * minimum is difficult (Picheny et al., 2012).
 * It is named after Howard H. Rosenbrock and was first published 1960 in "An automatic
 * method for finding the greatest or least value of a function. The Computer
 * Journal".
 *
 * Its optimal parameter = (1, ..., 1) and optimal function value = 0.
 *
 *         ⎛      ⎛                     ⎞²    ⎛          ⎞²  ⎞
 *     D-1 ⎜      ⎜                     ⎟     ⎜          ⎟   ⎟
 *      ∑  ⎜100 * ⎜ p(i + 1) - (p(i))²  ⎟   + ⎜ p(i) - 1 ⎟   ⎟
 *     i=1 ⎜      ⎜                     ⎟     ⎜          ⎟   ⎟
 *         ⎝      ⎝                     ⎠     ⎝          ⎠   ⎠
 */
class rosenbrock_function : public problem
{
public:
  /**
   * Initialises a rosenbrock function with `dimension` dimensions, lower bounds
   * of -2.048 and upper bounds of 2.048.
   */
  rosenbrock_function(const arma::uword dimension);

  virtual double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
