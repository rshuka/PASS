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
   * the optimisation configuration. (`inertia`,
   * `cognitive_acceleration`, `social_acceleration`,
   * `attraction_center`)
   */
  const pass::parallel_swarm_search &pss;

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
   * The personal best position that this particle ever found.
   */
  arma::vec personal_best_position;

  /**
   * The fitness value of `personal_best_position`.
   */
  double personal_best_fitness_value;

  /**
   * Initialises this particle with a random position within the search space of
   * `problem` and immediately evaluates it.
   */
  particle(const pass::parallel_swarm_search &pss,
           const pass::problem &problem);

  /**
   * Calculates new velocity and position values for this particle, and
   * evaluates the new position. If the new position yields a better fitness
   * value, updates `personal_best_position` and `personal_best_fitness_valu` and returns `true`.
   *
   * `local_best_agent` is used as the best previous position in the
   * neighbourhood.
   */
  bool update(const pass::particle &local_best_agent);
};
} // namespace pass
