#pragma once

#include "pass_bits/optimiser.hpp"

namespace pass {
/**
 * Implements the [Standard Particle Swarm Optimisation 2011]
 * (http://ieeexplore.ieee.org/xpls/icp.jsp?arnumber=6557848) algorithm.
 */
class parallel_swarm_search : public optimiser {
 public:
  /**
   * The velocity of all particles is initialized with a uniformly distributed
   * random value from the interval [-initial_velocity, initial_velocity].
   *
   * Is initialized to `0.5`.
   */
  double initial_velocity;

  /**
   * Controls the inertia of the particles:
   *
   *     new_velocity = old_velocity * maximal_acceleration + attraction_center
   *
   * Is initialized to `1 / (2 * log(2))`.
   */
  double maximal_acceleration;

  /**
   * Particles are pulled towards their personal best found parameter with a
   * random strength picked from the range [0, maximal_local_attraction].
   *
   * Is initialized to `0.5 + log(2)`.
   */
  double maximal_local_attraction;

  /**
   * Particles are pulled towards the global best found parameter with a random
   * strength picked from the range [0, maximal_global_attraction].
   *
   * Is initialized to `0.5 + log(2)`.
   */
  double maximal_global_attraction;

  /**
   * Used by the [agile restart mechanism](TODO: citation). The optimization is
   * aborted (so it can be restarted by the caller) when the globally best
   * found value better than the swarm average of historically best found times
   * `stagnationThreshold`.
   *
   * Is initialized to `0.9`. Should be set to a value between `0.9` for low-
   * dimensional problems and `0.7` for high-dimensional problems according to
   * the paper.
   */
  double stagnationThreshold;

  /**
   * The number of particles used during optimisation.
   *
   * Is initialized to `40`.
   */
  arma::uword population_size;

  parallel_swarm_search() noexcept;

  virtual optimise_result optimise(const pass::problem& problem);

 private:
  optimise_result optimise(const pass::problem& problem,
                           const std::chrono::nanoseconds maximal_duration);
};
}  // namespace pass
