#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * `messenger_full` is a 26-dimensional optimization problem issued by the ESA:
 * https://www.esa.int/gsp/ACT/projects/gtop/messenger_full.html
 */
class messenger_full : public problem
{
public:
  /**
   * Initializes the lower and upper bounds to the values listed on the ESA page.
   */
  messenger_full();

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
