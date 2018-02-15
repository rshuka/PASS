#pragma once

#include <armadillo>
#include <vector>
#include <cassert>
#include <algorithm> 
#include <cmath> 
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
   * bounds.
   */
  problem(const arma::uword dimension, const double lower_bound,
          const double upper_bound);

  /**
   * Initialises this problem with custom lower and upper bounds. The problem
   * dimension is derived from the bounds dimensions, which have to be equal.
   */
  problem(const arma::vec &lower_bounds, const arma::vec &upper_bounds);

  /**
   * Evaluates this problem at `parameter`, which must match the dimensions of
   * this problem.
   */
  virtual double evaluate(const arma::vec &parameter) const = 0;

  /**
   * Draws `count` uniformly distributed random vectors from range
   * [lower_bounds, upper_bounds].
   */
  arma::mat random_parameters(const arma::uword count) const;
};

}  // namespace pass
