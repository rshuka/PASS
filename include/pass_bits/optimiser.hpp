#pragma once

#include <chrono>

#include <armadillo>

#include "problem.hpp"

namespace pass {

/**
 * Created by `optimiser::optimise`. Stores information about the optimisation.
 */
struct optimise_result {
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
   * The total number of times `problem.evaluate` was called.
   */
  arma::uword evaluations;

  /**
   * The number of iterations executed by the optimiser.
   */
  arma::uword iterations;

  /**
   * Total time the optimiser took to find `objective_value`.
   */
  std::chrono::nanoseconds duration;

  optimise_result(const arma::uword dimension) noexcept;
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
   * The maximal number of evaluations before the optimiser stops.
   *
   * Initialised to the maximum value representable with `arma::uword`.
   */
  arma::uword maximal_evaluations;

  /**
   * The maximal elapsed time before the optimiser stops.
   *
   * Initialised to 1 second.
   */
  std::chrono::nanoseconds maximal_duration;

  optimiser() noexcept;

  /**
   * Optimises `problem`, starting at `initial_parameters`.
   *
   * `initial_parameters` must be a N×M matrix, with N = `problem.dimension()`
   * and M = number of parameters, which depends on the optimiser.
   */
  virtual optimise_result optimise(const pass::problem& problem,
                                   const arma::mat& initial_parameters) = 0;
};
}  // namespace pass
