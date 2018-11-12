#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * The Griewank function has many widespread local minima, which
 * are regularly distributed.
 *
 * Its optimal parameter = (0, ..., 0) and optimal function value = 0.
 *
 * \f[
 *   f(x_1 \cdots x_n) = 1 + \frac{1}{4000} \sum_{i=1}^n x^2_i - \prod_{i=1}^n cos(\frac{x_i}{\sqrt{i}})
 * \f]
 *
 *  \f[
 *   -600.00 \leq x_i \leq 6.00
 * \f]
 *   \f[
 *   \text{minimum at }f(0, \cdots, 0) = 0
 * \f]
 */

//
//     D                  D
//     Σ (p(i)² / 4000) - ∏ cos(p(i)/√(i)) + 1
//    i=1                i=1

class griewank_function : public problem
{
public:
  /**
   * Initialises a Griewank function with `dimension` dimensions, lower bounds of
   * -600.0 and upper bounds of 600.0.
   */
  explicit griewank_function(const arma::uword dimension);

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
