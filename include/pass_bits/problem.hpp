#pragma once

#include "pass_bits/config.hpp"
#include <algorithm>
#include <armadillo>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>

#if defined(SUPPORT_MPI)
#include <mpi.h>
#endif
namespace pass
{
/**
 * This class defines the interface for problems that can be approximated by a
 * `pass::optimiser`. Subclasses need only implement `evaluate`.
 */
class problem
{
public:
  /**
   * Lower bound constraints for each problem dimension. An optimiser will never
   * evaluate a point outside of these bounds.
   *
   * Must be element-wise less than `upper_bounds`.
   */
  const arma::vec lower_bounds;

  /**
   * Upper bound constraints for each problem dimension. An optimiser will never
   * evaluate a point outside of these bounds.
   *
   * Must be element-wise greater than `upper_bounds`.
   */
  const arma::vec upper_bounds;

  /**
   * Identify every problem with its own name
   */
  const std::string name;

  /**
   * Returns the element-wise maximum difference between `upper_bounds` and
   * `lower_bounds`.
   */
  arma::vec bounds_range() const noexcept;

  /**
   * Returns the problem dimension.
   */
  arma::uword dimension() const noexcept;

  /**
   * Initialises an `dimension`-dimensional problem with uniform lower and upper
   * bounds and a problem name.
   */
  problem(const arma::uword dimension, const double &lower_bound,
          const double &upper_bound, const std::string &name);

  /**
   * Initialises this problem with custom lower and upper bounds and a problem name.
   * The problem dimension is derived from the bounds dimensions, which have to be equal.
   */
  problem(const arma::vec &lower_bounds, const arma::vec &upper_bounds, const std::string &name);

  /**
   * Evaluates this problem at `agent`, which must match the dimensions of
   * this problem.
   */
  virtual double evaluate(const arma::vec &agent) const = 0;

  /**
   * Evaluates this problem at `agent`, which must be a normalized vector (all
   * values must be in range [0, 1]). `agent` is mapped to the problem
   * boundaries before evaluation.
   */
  double evaluate_normalised(const arma::vec &normalised_agent) const;

  /**
   * Draws `count` uniformly distributed random agents from range [0, 1], stored
   * column-wise.
   */
  arma::mat normalised_random_agents(const arma::uword count) const;

  /**
   * Generates `count` hammersley points with the same dimension as this, stored
   * column-wise.
   *
   * For a definition of hammersley points, see
   * http://www.cse.cuhk.edu.hk/~ttwong/papers/udpoint/udpoint.pdf,
   * equations 1 to 3.
   */
  arma::mat normalised_hammersley_agents(const arma::uword count) const;

  /**
   * Draws `count` distributed random agents from range [0, 1], stored
   * column-wise.
   * Agens are from type hammersley or random distributed (50%-50%)
   */
  arma::mat initialise_normalised_agents(const arma::uword count) const;
};

} // namespace pass
