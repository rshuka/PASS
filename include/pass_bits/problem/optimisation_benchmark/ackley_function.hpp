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
 * \f[
 *   f(x_0 \cdots x_n) = -20 exp(-0.2 \sqrt{\frac{1}{n} \sum_{i=1}^n x_i^2}) - exp(\frac{1}{n} \sum_{i=1}^n cos(2\pi x_i)) + 20 + e
 * \f]
 *
 *  \f[
 *   -32 \leq x_i \leq 32
 * \f]
 *   \f[
 *   \text{minimum at }f(0, \cdots, 0) = 0
 * \f]
 */

//           ⎛        ⎛ D      ⎞⎞       ⎛ D               ⎞
//           ⎜        ⎜ ∑ p(i)²⎟⎟       ⎜ ∑ cos(2π * p(i))⎟
// -20 * exp ⎜-0.2 * √⎜i=1     ⎟⎟ - exp ⎜i=1              ⎟ + 20 + exp(1)
//           ⎜        ⎜--------⎟⎟       ⎜-----------------⎟
//           ⎝        ⎝   D    ⎠⎠       ⎝        D        ⎠
//
//
class ackley_function : public problem
{
public:
  /**
   * Initialises a ackley function with `dimension` dimensions, lower bounds of
   * -32.768 and upper bounds of 32.768.
   */
  explicit ackley_function(const arma::uword dimension);

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
