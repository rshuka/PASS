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
 *                 D ⎛                        ⎞
 *  418.9829 * D - ∑ ⎜ p(i) * sin(√ ||p(i)||) ⎟
 *                i=1⎝                        ⎠
 */
class schwefel_function : public problem
{
public:
  /**
   * Initialises a schwefel function with `dimension` dimensions, lower bounds
   * of -500.0 and upper bounds of 500.0.
   */
  schwefel_function(const arma::uword dimension);

  virtual double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
