#include "pass_bits/optimiser/parallel_swarm_search/parallel_swarm_search.hpp"
#include "pass_bits/optimiser/parallel_swarm_search/particle.hpp"

// pass::random_number_generator(), pass::random_neighbour()
#include "pass_bits/helper/random.hpp"

// std::pow
#include <cmath>

// assert
#include <cassert>

// std::vector
#include <vector>

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

  auto start_time = std::chrono::steady_clock::now();

  // Initialisation of PSO
  pass::optimise_result result(problem.dimension());
  pass::particle* global_best;

  // Instantiate `population_size` particles.
  std::vector<particle> particles;
  particles.reserve(population_size);
  for (arma::uword n = 0; n < population_size; ++n) {
    particles.push_back({*this, problem});

    ++result.evaluations;
    result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start_time);

    if (particles[n].best_value <= result.objective_value) {
      result.parameter = particles[n].best_parameter;
      result.objective_value = particles[n].best_value;
      global_best = &particles[n];

      if (result.objective_value <= acceptable_objective_value) {
        result.solved = true;
        return result;
      }
    }

    if (result.evaluations >= maximal_evaluations) {
      return result;
    } else if (result.duration >= maximal_duration) {
      return result;
    }
  }
  ++result.iterations;

  while (result.duration < maximal_duration &&
         result.evaluations < maximal_evaluations &&
         result.objective_value > acceptable_objective_value) {
    pass::particle& particle = particles[result.evaluations % population_size];

    // If `particle.update()` returns `true`, it has found a new personal best
    // value.
    if (particle.update(*global_best)) {
      if (particle.best_value < result.objective_value) {
        result.parameter = particle.best_parameter;
        result.objective_value = particle.best_value;
        global_best = &particle;
      }
    }

    if (++result.evaluations == population_size) {
      ++result.iterations;
    }

    result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start_time);
  }

  result.solved = result.objective_value <= acceptable_objective_value;
  return result;
}
