#include "pass_bits/optimiser/parallel_swarm_search/parallel_swarm_search.hpp"
#include "pass_bits/optimiser/parallel_swarm_search/particle.hpp"
#include "pass_bits/helper/random.hpp"
#include "pass_bits/helper/stopwatch.hpp"
#include <cmath> // std::pow
#include <vector> // std::vector

pass::parallel_swarm_search::parallel_swarm_search() noexcept
    : optimiser(),
      initial_velocity(0.5),
      maximal_acceleration(1.0 / (2.0 * std::log(2.0))),
      maximal_local_attraction(0.5 + std::log(2.0)),
      maximal_global_attraction(maximal_local_attraction),
      population_size(40) {}

pass::optimise_result pass::parallel_swarm_search::optimise(const pass::problem& problem) 
{
  assert(initial_velocity >= 0.0);
  assert(maximal_acceleration >= 0.0);
  assert(maximal_local_attraction >= 0.0);
  assert(maximal_global_attraction >= 0.0);

  pass::stopwatch stopwatch;

  pass::optimise_result result(problem.dimension(), acceptable_objective_value,
                               population_size);

  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved()) 
  {
    pass::optimise_result sub_result =
        optimise(problem, maximal_duration - result.duration);
    if (sub_result.objective_value < result.objective_value) 
    {
      result.parameter = sub_result.parameter;
      result.objective_value = sub_result.objective_value;
    }
    result.iterations += sub_result.iterations;
    result.duration = stopwatch.get_elapsed();
  }
  return result;
}

pass::optimise_result pass::parallel_swarm_search::optimise(
    const pass::problem& problem,
    const std::chrono::nanoseconds maximal_duration) 
{
  // Note: `result.duration` isn't used by this function because the termination
  // criterion can be extracted from the stopwatch directly, and the caller has
  // its own timer anyways.
  pass::stopwatch stopwatch;

  pass::optimise_result result(problem.dimension(), acceptable_objective_value,
                               population_size);
  pass::particle* global_best = nullptr;

  // Instantiate `population_size` particles.
  std::vector<particle> particles;
  particles.reserve(population_size);
  for (arma::uword n = 0; n < population_size; ++n) 
  {
    particles.push_back({*this, problem});

    if (particles[n].best_value <= result.objective_value) 
    {
      result.parameter = particles[n].best_parameter;
      result.objective_value = particles[n].best_value;
      global_best = &particles[n];
    }
  }
  ++result.iterations;

  while (stopwatch.get_elapsed() < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved()) 
  {
    for (pass::particle& particle : particles) 
    {
      // If `particle.update()` returns `true`, it has found a new personal best
      // value.
      if (particle.update(*global_best)) 
      {
        if (particle.best_value < result.objective_value) 
        {
          result.parameter = particle.best_parameter;
          result.objective_value = particle.best_value;
          global_best = &particle;
        }
      }
    }

    ++result.iterations;

    // Restart criterion
    const double averagePerformance =
        std::accumulate(particles.begin(), particles.end(), 0.0,
                        [](const double sum, const pass::particle& particle) {
                          return sum + particle.best_value;
                        });
    if (result.objective_value >= stagnationThreshold * averagePerformance) 
    {
      break;
    }
  }

  return result;
}