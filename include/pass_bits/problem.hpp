#pragma once

#include <algorithm>
#include <armadillo>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>
#if defined(SUPPORT_OPENMP)
#include <omp.h>
#endif

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
  problem(const arma::uword dimension, const double lower_bound,
          const double upper_bound, const std::string name);

  /**
   * Initialises this problem with custom lower and upper bounds and a problem name.
   * The problem dimension is derived from the bounds dimensions, which have to be equal.
   */
  problem(const arma::vec &lower_bounds, const arma::vec &upper_bounds, const std::string name);

  /**
   * Evaluates this problem at `agent`, which must match the dimensions of
   * this problem.
   */
  virtual double evaluate(const arma::vec &agent) const = 0;

  /**
   * Draws `count` uniformly distributed random agents from range
   * [lower_bounds, upper_bounds], stored column-wise.
   */
  arma::mat random_agents(const arma::uword count) const;

  /**
   * Generates `count` hammersley points with the same dimension as this,
   * mapped to the problem bounds of this, stored column-wise.
   *
   * For a definition of hammersley points, see
   * http://www.cse.cuhk.edu.hk/~ttwong/papers/udpoint/udpoint.pdf,
   * equations 1 to 3.
   */
  arma::mat hammersley_agents(const arma::uword count) const;
};

} // namespace pass
