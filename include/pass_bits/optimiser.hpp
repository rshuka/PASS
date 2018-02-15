#pragma once

#include <chrono> // std:: chrono
#include <armadillo> // std::arma::vec
#include <cassert> // assert
#include "pass_bits/problem.hpp"

#if defined(SUPPORT_OPENMP)
#include <omp.h> 
#endif

#if defined(SUPPORT_MPI)
#include <mpi.h> 
#endif

namespace pass 
{
/**
 * Created by `optimiser::optimise`. Stores information about the optimisation.
 */
struct optimise_result 
{
  /**
   * The parameter that was evaluated to `objective_value`. Is initialised to
   * NaN in all dimensions.
   */
  arma::vec parameter;

  /**
   * The best found objective value. Is initialised to infinity.
   */
  double objective_value;

  /**
   * The acceptable objective value of the optimiser that created this result.
   */
  const double acceptable_objective_value;

  /**
   * The number of iterations executed by the optimiser.
   */
  arma::uword iterations;

  /**
   * Total time the optimiser took to find `objective_value`.
   */
  std::chrono::nanoseconds duration;

  /**
   * The total number of times `problem.evaluate` was called per [iteration].
   */
  const arma::uword evaluations_per_iteration;

  /**
   * Initializes [evaluations_per_iteration] to `1`.
   */
  optimise_result(const arma::uword dimension,
                  const double acceptable_objective_value) noexcept;

  optimise_result(const arma::uword dimension,
                  const double acceptable_objective_value,
                  const arma::uword evaluations_per_iteration) noexcept;

  /**
   * Returns `true` if the problem was solved.
   */
  bool solved() const;
};

/**
 * Interface for algorithms that approximate a solution to a `pass::problem`.
 * Subclasses need only implement `optimise`.
 */
class optimiser {
 public:
  /**
   * The optimiser will stop early if it finds a parameter that evaluates to a
   * value less than or equal to `acceptable_objective_value`.
   *
   * Initialised to -infinity.
   */
  double acceptable_objective_value;

  /**
   * The maximal number of iterations before the optimiser stops.
   *
   * Initialised to the maximum value representable with `arma::uword`.
   */
  arma::uword maximal_iterations;

  /**
   * The maximal elapsed time before the optimiser stops.
   *
   * Initialised to 1 minute.
   */
  std::chrono::nanoseconds maximal_duration;

  optimiser() noexcept;

  /**
   * Optimises `problem`, storing the result and performance characteristics of
   * the optimisation in the returned `optimise_result`.
   */
  virtual optimise_result optimise(const pass::problem& problem) = 0;
};
}  // namespace pass
