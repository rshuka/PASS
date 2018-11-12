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
 *
 * \f[
 *   f(x_1 \cdots x_n) = \frac{1}{2} \sum_{i=1}^n (x_i^4 - 16x_i^2 + 5x_i)
 * \f]
 *
 *  \f[
 *   -5.00 \leq x_i \leq 5.00
 * \f]
 *   \f[
 *   \text{minimum at }f(-2.903534, \cdots, -2.903534) = -39.16599 \cdot n
 * \f]
 */

//
//          D ⎛                             ⎞
//   0.5  * ∑ ⎜p(i)⁴ - 16 * p(i)² + 5 * p(i)⎟
//         i=1⎝                             ⎠
//
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
