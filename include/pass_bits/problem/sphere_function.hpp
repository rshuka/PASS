#pragma once

#include "../problem.hpp"

namespace pass {

/**
 * The sphere function is a common *toy* problem with a very small computational
 * cost, used for testing and benchmarking algorithms.
 *
 *      n
 *      Σ (p(i)^2)
 *     i=1
 */
class sphere_function : public problem {
 public:
  virtual double evaluate(const arma::vec& parameter) const override;

  /**
   * Initialises a sphere function with `dimension` dimensions, lower bounds of
   * -5.12 and upper bounds of 5.12.
   */
  sphere_function(const arma::uword dimension);
};
}  // namespace pass
