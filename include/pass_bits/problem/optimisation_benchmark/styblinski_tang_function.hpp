#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * The styblinski function is a common *toy* problem with a very small computational
 * cost, used for testing and benchmarking algorithms. It is a multimodal function.
 *
 * Its optimal parameter = (-2.903534, ..., -2.903534) and optimal function value = -39.16599 * Dimension.
 *
 *           D ⎛                             ⎞
 *    0.5  * ∑ ⎜p(i)⁴ - 16 * p(i)² + 5 * p(i)⎟
 *          i=1⎝                             ⎠
 */
class styblinski_tang_function : public problem
{
public:
  /**
   * Initialises a styblinski function with `dimension` dimensions, lower bounds of
   * -5.0 and upper bounds of 5.0.
   */
  explicit styblinski_tang_function(const arma::uword dimension);

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
