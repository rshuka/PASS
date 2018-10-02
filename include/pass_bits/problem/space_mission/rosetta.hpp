#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * `rosetta` is a 22-dimensional optimization problem issued by the ESA:
 * https://www.esa.int/gsp/ACT/projects/gtop/rosetta.html
 */
class rosetta : public problem
{
public:
  /**
   * Initializes the lower and upper bounds to the values listed on the ESA page.
   */
  rosetta();

  virtual double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
