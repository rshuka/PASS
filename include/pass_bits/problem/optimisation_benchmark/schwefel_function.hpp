#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * The Schwefel function is complex, with many local minima.
 *
 * Its optimal parameter = (420.9687, ..., 420.9687) and optimal function value = 0.
 * In this implementation for the optimal parameter the function value = 5.42627e-10.
 *
 * \f[
 *   f(x_1 \cdots x_n) = \sum_{i=1}^n (-x_i sin(\sqrt{|x_i|})) + \alpha \cdot n
 * \f]
 *
 * \f[
 * \alpha = 418.982887
 * \f]
 *
 * \f[
 *   500.00 \leq x_i \leq 500.00
 * \f]
 *
 * \f[
 *   \text{minimum at }f(420.968746, 420.968746, \cdots, 420.968746) = 0
 * \f]
 */

//
//                         D ⎛                        ⎞
// 418.9828872724338 * D - ∑ ⎜ p(i) * sin(√ ||p(i)||) ⎟
//                        i=1⎝                        ⎠
//
class schwefel_function : public problem
{
public:
  /**
   * Initialises a schwefel function with `dimension` dimensions, lower bounds
   * of -500.0 and upper bounds of 500.0.
   */
  explicit schwefel_function(const arma::uword dimension);

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
