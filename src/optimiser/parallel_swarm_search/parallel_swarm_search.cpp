#include "pass_bits/optimiser/parallel_swarm_search/parallel_swarm_search.hpp"
#include "pass_bits/optimiser/parallel_swarm_search/particle.hpp"
#include "pass_bits/helper/random.hpp"
#include "pass_bits/helper/stopwatch.hpp"
#include <cmath>  // std::pow
#include <vector> // std::vector

pass::parallel_swarm_search::parallel_swarm_search() noexcept
    : optimiser(),
      swarm_size(40),
      inertia(1.0 / (2.0 * std::log(2.0))),
      cognitive_acceleration(0.5 + std::log(2.0)),
      social_acceleration(cognitive_acceleration),
      neighbourhood_probability(1.0 -
                                std::pow(1.0 - 1.0 / static_cast<double>(swarm_size), 3.0)),
      stagnation_threshold(0.9) {}

pass::optimise_result pass::parallel_swarm_search::optimise(
    const pass::problem &problem)
{
  assert(inertia >= 0.0);
  assert(cognitive_acceleration >= 0.0);
  assert(social_acceleration >= 0.0);
  assert(neighbourhood_probability > 0.0 && neighbourhood_probability <= 1.0);

  pass::stopwatch stopwatch;
  stopwatch.start();

  pass::optimise_result result(problem.dimension(), acceptable_fitness_value);

  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved())
  {
    pass::optimise_result sub_result =
        optimise(problem, maximal_duration - result.duration);

    if (sub_result.fitness_value < result.fitness_value)
    {
      result.agent = sub_result.agent;
      result.fitness_value = sub_result.fitness_value;
    }

    result.iterations += sub_result.iterations;
    result.duration = stopwatch.get_elapsed();
  }
  return result;
}

pass::optimise_result pass::parallel_swarm_search::optimise(
    const pass::problem &problem,
    const std::chrono::nanoseconds maximal_duration)
{
  // Note: `result.duration` isn't used by this function because the termination
  // criterion can be extracted from the stopwatch directly, and the caller has
  // its own timer anyways.
  pass::stopwatch stopwatch;
  stopwatch.start();

  pass::optimise_result result(problem.dimension(), acceptable_fitness_value);

  // Instantiate `swarm size` particles.
  std::vector<particle> particles;
  particles.reserve(swarm_size);

  for (arma::uword n = 0; n < swarm_size; ++n)
  {
    particles.push_back({*this, problem});

    if (particles[n].personal_best_fitness_value <= result.fitness_value)
    {
      result.agent = particles[n].personal_best_position;
      result.fitness_value = particles[n].personal_best_fitness_value;
    }
  }
  ++result.iterations;

  arma::umat topology(swarm_size, swarm_size);
  bool randomize_topology = true;

  while (stopwatch.get_elapsed() < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved())
  {
    if (randomize_topology)
    {
      topology = (arma::mat(swarm_size, swarm_size,
                            arma::fill::randu) < neighbourhood_probability);
      // When searching for the best neighbour, we begin with the particles
      // personal best value; We don't need to visit it twice.
      topology.diag().fill(0);
    }
    randomize_topology = true;

    for (arma::uword n = 0; n < swarm_size; n++)
    {
      pass::particle &particle = particles[n];

      pass::particle *local_best_agent = &particle;

      for (arma::uword i = 0; i < swarm_size; i++)
      {
        if (topology(n, i) && particles[i].personal_best_fitness_value < local_best_agent->personal_best_fitness_value)
        {
          local_best_agent = &particles[i];
        }
      }

      // If `particle.update()` returns `true`, it has found a new personal best
      // value.
      if (particle.update(*local_best_agent))
      {
        if (particle.personal_best_fitness_value < result.fitness_value)
        {
          result.agent = particle.personal_best_position;
          result.fitness_value = particle.personal_best_fitness_value;
          randomize_topology = false;
        }
      }
    }

    ++result.iterations;

    continue;
    // Restart criterion
    const double average_performance =
        std::accumulate(particles.begin(), particles.end(), 0.0,
                        [](const double sum, const pass::particle &particle) {
                          return sum + particle.personal_best_fitness_value;
                        });
    if (result.fitness_value >= stagnation_threshold * average_performance)
    {
      break;
    }
  }

  return result;
}
