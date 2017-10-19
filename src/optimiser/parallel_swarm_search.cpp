#include "pass_bits/optimiser/parallel_swarm_search.hpp"

// pass::random_number_generator(), pass::random_neighbour()
#include "pass_bits/helper/random.hpp"

// std::pow
#include <cmath>

// assert
#include <cassert>

pass::parallel_swarm_search::parallel_swarm_search() noexcept
    : optimiser(),
      initial_velocity(0.5),
      maximal_acceleration(1.0 / (2.0 * std::log(2.0))),
      maximal_local_attraction(0.5 + std::log(2.0)),
      maximal_global_attraction(maximal_local_attraction),
      population_size(40) {}

pass::optimise_result pass::parallel_swarm_search::optimise(
    const pass::problem& problem) {
  assert(initial_velocity >= 0.0);
  assert(maximal_acceleration >= 0.0);
  assert(maximal_local_attraction >= 0.0);
  assert(maximal_global_attraction >= 0.0);

  return 0;
}
