#pragma once

#include "pass_bits/problem.hpp"

namespace pass
{

/**
 * This helper class reduces the problem boundaries of a wrapped problem to
 * 1/nth of its normal search space, where `n` can be specified as a constructor
 * argument.
 */
class search_space_constraint : public problem
{
public:
  /**
   * The internal problem to which `evaluate` calls are forwarded.
   *
   * CAUTION: This variable is a reference type! Do not free or reuse the
   * memory of `wrapped_object` until `search_space_constraint` is destroyed!
   */
  const pass::problem &wrapped_problem;

  /**
   * Initializes the problem boundaries. All boundaries are copied from
   * `wrapped_problem`, except for the dimension 1.
   *
   * The dimension 1 boundary is modified to have a range of
   * `1 / total_segments * original bounds range`, and and offset of
   * `calculated range * segment`.
   *
   * Both parameters `segment` and `total_segments` are 1-based: If you want to
   * divide the search space into four non-overlapping sections, pass the values
   * 1, 2, 3, 4 as `segment` parameters, and `4` as `total_segments` parameter.
   */
  search_space_constraint(const pass::problem &wrapped_problem,
                          arma::uword segment, arma::uword total_segments);

  double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
