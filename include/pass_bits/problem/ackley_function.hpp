#pragma once

#include "../problem.hpp"

namespace pass {

/**
 * The Ackley function is a common *toy* problem with a very small computational
 * cost, used for testing and benchmarking algorithms. It is named after David
 * H. Ackley and was first published 1987 in *A Connectionist Machine for
 * Genetic Hillclimbing. Kluwer Academic Publishers*. Its optimal parameter =
 * (0, ..., 0) and optimal function value = 0.
 *
 *          ⎛        ⎛ N      ⎞⎞      ⎛ N               ⎞
 *          ⎜        ⎜ ∑ p(i)²⎟⎟      ⎜ ∑ cos(2π * p(i))⎟
 * -20 * exp⎜-0.2 * √⎜i=1     ⎟⎟ - exp⎜i=1              ⎟ + 20 + exp(1)
 *          ⎜        ⎜--------⎟⎟      ⎜-----------------⎟
 *          ⎝        ⎝   N    ⎠⎠      ⎝        N        ⎠
 */
class ackley_function : public problem {
 public:
  /**
   * Initialises a sphere function with `dimension` dimensions, lower bounds of
   * -32.768 and upper bounds of 32.768.
   */
  ackley_function(const arma::uword dimension);

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
