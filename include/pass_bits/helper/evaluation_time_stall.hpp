#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{

/**
 * This helper class wraps another `problem` object and artifically increases
 * the calculation time. Every time `evaluation_time_stall::evaluate()` is
 * called, the call is forwarded to `wrapped_problem` and repeated
 * `repetitions` times.
 */
class evaluation_time_stall : public problem
{
public:
  /**
   * The internal problem to which `evaluate` calls are forwarded.
   *
   * CAUTION: This variable is a reference type! Do not free or reuse the
   * memory of `wrapped_object` until `evaluation_time_stall` is destroyed!
   */
  const pass::problem &wrapped_problem;

  /**
   * The number of times every `evaluate()` call is repeated. Must be 1 or
   * greater.
   */
  arma::uword repetitions;

  /**
   * Initializes this object with the same bounds as `wrapped_problem` and
   * `repetitions` set to 1.
   */
  explicit evaluation_time_stall(const pass::problem &wrapped_problem);

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
