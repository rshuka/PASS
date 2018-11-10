#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * Ackley’s is a widely used multimodal test function. It is named after David
 * H. Ackley and was first published 1987 in "A Connectionist Machine for
 * Genetic Hillclimbing. Kluwer Academic Publishers".
 *
 * Its optimal parameter = (0, ..., 0) and optimal function value = 0.
 *
 *          ⎛        ⎛ D      ⎞⎞      ⎛ D               ⎞
 *          ⎜        ⎜ ∑ p(i)²⎟⎟      ⎜ ∑ cos(2π * p(i))⎟
 * -20 * exp⎜-0.2 * √⎜i=1     ⎟⎟ - exp⎜i=1              ⎟ + 20 + exp(1)
 *          ⎜        ⎜--------⎟⎟      ⎜-----------------⎟
 *          ⎝        ⎝   D    ⎠⎠      ⎝        D        ⎠
 */
class ackley_function : public problem
{
public:
  /**
   * Initialises a ackley function with `dimension` dimensions, lower bounds of
   * -32.768 and upper bounds of 32.768.
   */
  explicit ackley_function(const arma::uword dimension);

  virtual double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
