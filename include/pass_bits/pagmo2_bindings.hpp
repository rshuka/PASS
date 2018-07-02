#pragma once

#include <pagmo/pagmo.hpp>

#include "pass_bits/optimiser.hpp"
#include "pass_bits/problem.hpp"

namespace pass
{
namespace pagmo2
{

// ----------------------------------------------
// problem adapter
// ----------------------------------------------

/**
 * This class is a pagmo2 compatible UDP (user defined problem). pagmo requires
 * at least two methods of its problem type: `fitness()`, which acts the same as
 * `problem::evaluate()` in PASS, and `get_bounds()`, which returns the lower
 * and upper bounds as a tuple. This class forwards all calls to these methods
 * to `wrapped_problem`, making arbitrary PASS problems compatible to pagmo.
 */
class problem_adapter
{
public:
  const pass::problem *wrapped_problem;

  problem_adapter();
  problem_adapter(const pass::problem &wrapped_problem);

  /**
   * Returns the fitness value of agent. pagmo2 supports equality and inequality
   * constraints, which would also be returned from this method call. Because
   * PASS doesn't contain any problems that use these constraints, the result
   * will always have length 1, which is the fitness value.
   */
  pagmo::vector_double fitness(const pagmo::vector_double &agent) const;

  /**
   * Returns a (lower_bounds, upper_bounds) tuple.
   */
  std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const;
};

// ----------------------------------------------
// algorithm adapter
// ----------------------------------------------

/**
 * Wrapper an optimiser from the pagmo2 library:
 * https://esa.github.io/pagmo2/docs/algorithm_list.html
 *
 * This is the super class for other pagmo algorithm adapters. Subclasses should
 * only implement `get_algorithm()`.
 */
class algorithm_adapter : public optimiser
{
public:
  /**
   * Must be at least 5 or pagmo will throw an exception. Is initialized to 0
   * because you probably want to use a higher number than that.
   */
  arma::uword population_size;

  /**
   * Initialises the optimiser with its name.
   */
  algorithm_adapter(const std::string &name) noexcept;

  /**
   * Obtains an algorithm object from `get_algorithm()`.
   * Repeatedly lets the pagmo algorithm evolve the population a single
   * iteration. Increasing the number of evolutions per iteration would probably
   * speed up the algorithm, but it would also increase the delay between
   * reaching a termination criteria and actually stopping.
   */
  virtual optimise_result optimise(const pass::problem &problem);

protected:
  /**
   * Constructs a pagmo algorithm with the default values, except that this
   * method circumvents the pagmo termination criteria because we have our own.
   */
  virtual pagmo::algorithm get_algorithm() const = 0;
};

// ----------------------------------------------
// CMAES
// ----------------------------------------------

/**
 * Wrapper for the CMAES optimiser from the pagmo2 library:
 * https://esa.github.io/pagmo2/docs/cpp/algorithms/cmaes.html
 *
 * Enforces evaluations to stay within the problem boundaries (even though the
 * pagmo documentation states it would make algorithm less efficient) because
 * that's what the PASS problem guarantees.
 */
class cmaes : public algorithm_adapter
{
public:
  /**
   * Initialises the optimiser with its name.
   */
  cmaes() noexcept;

protected:
  virtual pagmo::algorithm get_algorithm() const;
};

// ----------------------------------------------
// differential evolution
// ----------------------------------------------

/**
 * Wrapper for the differential evolution optimiser from the pagmo2 library:
 * https://esa.github.io/pagmo2/docs/cpp/algorithms/cmaes.html
 */
class differential_evolution : public algorithm_adapter
{
public:
  /**
   * Initialises the optimiser with its name.
   */
  differential_evolution() noexcept;

protected:
  virtual pagmo::algorithm get_algorithm() const;
};
} // namespace pagmo2
} // namespace pass
