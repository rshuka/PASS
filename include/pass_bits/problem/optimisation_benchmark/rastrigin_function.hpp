#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * Rastrigin's function is based on the function of De Jong with the addition
 * of cosine modulation to produce many local minima. Thus, the test function is highly
 * multimodal. However, the location of the minima are regularly distributed.
 * It is named after L. A. Rastrigin and was first published 1974 in "Systems of Extremal
 * Control" as 2-dimensional problem and 1991 generalised for n dimensions by H.
 * Mühlenbein, D. Schomisch and J. Born in "The Parallel Genetic Algorithm as
 * Function Optimiser. Parallel Computing".
 *
 * Its optimal parameter = (0, ..., 0) and optimal function value = 0.
 *
 *              D ⎛                            ⎞
 *     10 * D + ∑ ⎜p(i)² - 10 * cos(2π * p(i))⎟
 *             i=1⎝                            ⎠
 */
class rastrigin_function : public problem
{
public:
  /**
   * Initialises a rastrigin function with `dimension` dimensions, lower bounds
   * of -5.12 and upper bounds of 5.12.
   */
  explicit rastrigin_function(const arma::uword dimension);

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
