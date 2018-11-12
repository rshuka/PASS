#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * So called first function of De Jong’s is one of the simplest
 * test benchmark. Function is continuous, convex and unimodal.
 *
 * Its optimal parameter = (0, ..., 0) and optimal function value = 0.
 *
 * \f[
 *   f(x_1 \cdots x_n) = \sum_{i=1}^n x_i^2
 * \f]
 *
 *  \f[
 *   -5.12 \leq x_i \leq 5.12
 * \f]
 *   \f[
 *   \text{minimum at }f(0, \cdots, 0) = 0
 * \f]
 */

//
//     D
//     Σ (p(i)²)
//    i=1

class de_jong_function : public problem
{
public:
  /**
   * Initialises a De Jong's function with `dimension` dimensions, lower bounds of
   * -5.12 and upper bounds of 5.12.
   */
  explicit de_jong_function(const arma::uword dimension);

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
