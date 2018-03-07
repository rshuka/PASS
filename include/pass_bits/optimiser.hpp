#pragma once

#include "pass_bits/problem.hpp"
#include "pass_bits/helper/stopwatch.hpp"
#include <chrono>    // std::chrono
#include <armadillo> // std::arma::vec, arma::uword
#include <cassert>   // assert
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
 * Created by `optimiser::optimise`. Stores information about the optimisation.
 */
struct optimise_result
{
  /**
   * The best found parameter that was evaluated to `fitness_value`. Is initialised to
   * NaN in all dimensions.
   */
  arma::vec agent;

  /**
   * The best found objective value. Is initialised to + infinity.
   */
  double fitness_value;

  /**
   * The acceptable objective (good enough) value of the optimiser.
   */
  const double acceptable_fitness_value;

  /**
   * The number of iterations executed by the optimiser.
   */
  arma::uword iterations;

  /**
   * Total time in nanoseconds (10^-9) the optimiser took to find `objective_value`.
   */
  std::chrono::nanoseconds duration;

  optimise_result(const arma::uword dimension,
                  const double acceptable_fitness_value) noexcept;

  /**
   * Returns `true` if the problem was solved.
   */
  bool solved() const;
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
   * Initialised to 1 minute.
   */
  std::chrono::microseconds maximal_duration;

  optimiser() noexcept;

  /**
   * Optimises `problem`, storing the result and performance characteristics of
   * the optimisation in the returned `optimise_result`.
   */
  virtual optimise_result optimise(const pass::problem &problem) = 0;
};

} // namespace pass

/**
 * Overloading the `<<` operator enables the optimise_result to be printed to
 * the console.
 *
 * See: https://docs.microsoft.com/en-us/cpp/standard-library/overloading-the-output-operator-for-your-own-classes
 */
std::ostream& operator<<(std::ostream& os, const pass::optimise_result& result);
