#pragma once

#include "pass_bits/optimiser.hpp"

namespace pass
{
/**
 * Implements "Parallel Swarm Search" - a new Algorithm basen on
 * the Standard Particle Swarm Optimisation 2011 algorithm
 * (http://ieeexplore.ieee.org/xpls/icp.jsp?arnumber=6557848)
 *
 * The initialisation variables can be found in
 * (http://clerc.maurice.free.fr/pso/SPSO_descriptions.pdf)
 *
 */
class parallel_swarm_search : public optimiser
{
public:
  /**
   * The number of particles used during optimisation.
   *
   * Is initialized to `40`.
   */
  arma::uword swarm_size;

  /**
   * Controls the inertia of the particles:
   *
   * Is initialized to `1 / (2 * log(2))`.
   */
  double inertia;

  /**
   * Particles are pulled towards their personal best found parameter with a
   * strength
   *
   * Is initialized to `0.5 + log(2)`.
   */
  double cognitive_acceleration;

  /**
   * Particles are pulled towards the global best found parameter with a
   * strength
   *
   * Is initialized to `0.5 + log(2)`.
   */
  double social_acceleration;

  /**
   * Probability of a particle to be in the neighbourhood of another particle.
   * Must be in range `[0, 1]`.
   *
   * Is initialized to `1 - (1 - 1/population_size)^3`.
   * Based on Clerc description (Method 2)
   * http://clerc.maurice.free.fr/pso/random_topology.pdf
   */
  double neighbourhood_probability;

  /**
   * Denotes the migration invervall for the MPI Communication
   *
   * Is initialized to 0 i.e. MPI after every iteration
   */
#if defined(SUPPORT_MPI)
  arma::uword migration_stall;
#endif

  parallel_swarm_search() noexcept;

  virtual optimise_result optimise(const pass::problem &problem);

private:
  /**
   * The number of threads if openMP is enabled
   *
   * Is initialized to maximum available threads
   */
#if defined(SUPPORT_OPENMP)
  int number_threads;
#endif
};
} // namespace pass
