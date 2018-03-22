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
 *      D                  D
 *      Σ (p(i)² / 4000) - ∏ cos(p(i)/√(i)) + 1
 *     i=1                i=1
 */
class griewank_function : public problem
{
public:
  /**
   * Initialises a Griewank function with `dimension` dimensions, lower bounds of
   * -600.0 and upper bounds of 600.0.
   */
  griewank_function(const arma::uword dimension);

  virtual double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
