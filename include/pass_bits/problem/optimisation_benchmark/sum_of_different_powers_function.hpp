#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * The Sum of Different Powers function is unimodal.
 *
 * Its optimal parameter = (0, ..., 0) and optimal function value = 0.
 *
 * \f[
 *   f(x_1 \cdots x_n) = \sum_{i=1}^n |x_i|^{(i+1)}
 * \f]
 *
 *  \f[
 *   -1.00 \leq x_i \leq 1.00
 * \f]
 *   \f[
 *   \text{minimum at }f(0, \cdots, 0) = 0
 * \f]
 */

//
//     D  ⎛                 ⎞
//     ∑  ⎜||p(i)||^(i + 1) ⎟
//    i=1 ⎝                 ⎠
//
class sum_of_different_powers_function : public problem
{
public:
  /**
   * Initialises a different powers function with `dimension` dimensions, lower
   * bounds of -1.0 and upper bounds of 1.0.
   */
  explicit sum_of_different_powers_function(const arma::uword dimension);

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
