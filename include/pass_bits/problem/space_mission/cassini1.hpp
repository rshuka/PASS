#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * `cassini1` is a 6-dimensional optimization problem issued by the ESA:
 * https://www.esa.int/gsp/ACT/projects/gtop/cassini1.html
 */
class cassini1 : public problem
{
public:
  /**
   * Initializes the lower and upper bounds to the values listed on the ESA page.
   */
  cassini1();

  virtual double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
