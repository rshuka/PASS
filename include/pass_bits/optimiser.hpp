#pragma once

#include "pass_bits/problem.hpp"
#include "pass_bits/helper/stopwatch.hpp"
#include "pass_bits/config.hpp"
#include <stdexcept> // throw error
#include <chrono>    // std::chrono
#include <armadillo> // std::arma::vec, arma::uword
#include <cassert>   // assert
#include <string>    // name

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
   * The best found parameter that was evaluated to `fitness_value`. Is initialised to
   * NaN in all dimensions. This vector is mapped to the range [0, 1] in all dimensions.
   * To access the agent in problem search space coordinates, use
   * `optimise_result::agent()`.
   */
  arma::vec normalised_agent;

  /**
   * The best found objective value. Is initialised to + infinity.
   */
  double fitness_value;

  /**
   * The acceptable objective (good enough) value of the optimiser.
   */
  const double acceptable_fitness_value;

  /**
   * The problem that was optimised.
   */
  const pass::problem &problem;

  /**
   * The number of iterations executed by the optimiser.
   */
  arma::uword iterations;

  /**
   * The total number of times `problem.evaluate` was called.
   */
  arma::uword evaluations;

  /**
   * Total time in microseconds (10^-6) the optimiser took to find `objective_value`.
   */
  std::chrono::microseconds duration;

  optimise_result(const pass::problem &problem,
                  const double acceptable_fitness_value) noexcept;

  /**
   * Returns `true` if the problem was solved.
   */
  bool solved() const;

  /**
   * Returns `normalised_agent`, but mapped to the search space of `problem`.
   */
  arma::vec agent() const;
};

/**
 * Interface for algorithms that approximate a solution to a `pass::problem`.
 * Subclasses need only implement `optimise`.
 */
class optimiser
{
public:
  /**
   * The optimiser will stop early if it finds a agent that evaluates to a
   * value less than or equal to `acceptable_fitness_value`.
   *
   * Initialised to -infinity.
   */
  double acceptable_fitness_value;

  /**
   * The maximal number of iterations before the optimiser stops.
   *
   * Initialised to the maximum value representable with `arma::uword`.
   */
  arma::uword maximal_iterations;

  /**
   * The maximal elapsed time before the optimiser stops.
   *
   * Initialised maximum time std::chrono::system_clock::duration::max().
   */
  std::chrono::microseconds maximal_duration;

  /**
   * Identify every optimiser with its own name
   */
  const std::string name;

  /**
   * Initialises an optimiser with its own name
   */
  optimiser(const std::string name);

  /**
   * Optimises `problem`, storing the result and performance characteristics of
   * the optimisation in the returned `optimise_result`.
   */
  virtual optimise_result optimise(const pass::problem &problem) = 0;
};

} // namespace pass
