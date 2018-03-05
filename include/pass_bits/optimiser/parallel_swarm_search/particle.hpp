#pragma once

#include "parallel_swarm_search.hpp"

namespace pass
{
/**
 * Represents a single particle during a PSO run. Stores all properties specific
 * to the particle, and computes its own next position and velocity when
 * `update()` is called.
 */
class particle
{
public:
  /**
   * Reference to the PSO instance that created this particle. Is used to access
   * the optimisation configuration. (`initial_velocity`,
   * `maximal_acceleration`, `maximal_local_attraction`,
   * `maximal_global_attraction`)
   */
  const pass::parallel_swarm_search &pso;

  /**
   * Reference to the problem to optimise. Is used to access the search space
   * boundaries, and to call `evaluate()` on.
   */
  const pass::problem &problem;

  /**
   * The current position. Can't be outside the problems search space.
   */
  arma::vec position;

  /**
   * The current velocity.
   */
  arma::vec velocity;

  /**
   * The personal best parameter that this particle ever found.
   */
  arma::vec best_agent;

  /**
   * The objective value of `best_parameter`.
   */
  double best_value;

  /**
   * Initialises this particle with a random position within the search space of
   * `problem` and immediately evaluates it.
   */
  particle(const pass::parallel_swarm_search &pso,
           const pass::problem &problem);

  /**
   * Calculates new velocity and position values for this particle, and
   * evaluates the new position. If the new position yields a better objective
   * value, updates `best_parameter` and `best_value` and returns `true`.
   *
   * `best_neighbour` is used as the best previous position in the
   * neighbourhood.
   */
  bool update(const pass::particle &best_neighbour);
};
} // namespace pass
