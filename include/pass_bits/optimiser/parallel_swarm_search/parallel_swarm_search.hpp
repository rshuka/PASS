#pragma once

// Armadillo
#include <armadillo>

// OpenMP
#if defined(SUPPORT_OPENMP)
#include <omp.h>
#endif

// PASS
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
   * The number of particles used during optimisation.
   *
   * Is initialized to `40`.
   */
  arma::uword population_size;

  parallel_swarm_search() noexcept;

  virtual optimise_result optimise(const pass::problem& problem);
};
}  // namespace pass