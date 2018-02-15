#pragma once

#include "pass_bits/optimiser.hpp"

namespace pass 
{
/**
 * Implements the [Standard Particle Swarm Optimisation 2011]
 * (http://ieeexplore.ieee.org/xpls/icp.jsp?arnumber=6557848) algorithm.
 */
class particle_swarm_optimisation : public optimiser {
 public:
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
   * The number of particles used during optimisation.
   *
   * Is initialized to `40`.
   */
  arma::uword population_size;

  /**
   * Probability of a particle to be in the neighbourhood of another particle.
   * Must be in range `[0, 1]`.
   *
   * Is initialized to `(1 - 1/population_size)^3`.
   */
  double neighbourhood_probability;

  particle_swarm_optimisation() noexcept;

  virtual optimise_result optimise(const pass::problem& problem);
};
}  // namespace pass
