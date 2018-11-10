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
 *      D
 *      Σ (p(i)²)
 *     i=1
 */
class de_jong_function : public problem
{
public:
  /**
   * Initialises a De Jong's function with `dimension` dimensions, lower bounds of
   * -5.12 and upper bounds of 5.12.
   */
  explicit de_jong_function(const arma::uword dimension);

  virtual double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
