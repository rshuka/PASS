#include "pass_bits/optimiser/particle_swarm_optimisation.hpp"

// pass::random_number_generator(), pass::random_neighbour()
#include "pass_bits/helper/random.hpp"
#include "pass_bits/helper/stopwatch.hpp"

// std::pow
#include <cmath>

// assert
#include <cassert>

pass::particle_swarm_optimisation::particle_swarm_optimisation() noexcept
    : optimiser(),
      maximal_acceleration(1.0 / (2.0 * std::log(2.0))),
      maximal_local_attraction(0.5 + std::log(2.0)),
      maximal_global_attraction(maximal_local_attraction),
      population_size(40),
      neighbourhood_probability(
          std::pow(1.0 - 1.0 / static_cast<double>(population_size), 3.0)) {}

pass::optimise_result pass::particle_swarm_optimisation::optimise(
    const pass::problem& problem) {
  assert(maximal_acceleration >= 0.0);
  assert(maximal_local_attraction >= 0.0);
  assert(maximal_global_attraction >= 0.0);
  assert(neighbourhood_probability > 0.0 && neighbourhood_probability <= 1.0);

  pass::stopwatch stopwatch;

  // Initialisation of PSO
  pass::optimise_result result(problem.dimension(), acceptable_objective_value,
                               population_size);

  // Particle data, stored column-wise.
  arma::mat positions = problem.random_parameters(population_size);

  arma::mat velocities(problem.dimension(), population_size, arma::fill::randu);

  // Place the velocities within the boundaries
  velocities.each_col() %= problem.bounds_range();
  velocities.each_col() += problem.lower_bounds;
  velocities -= positions;

  arma::mat best_found_parameters = positions;
  arma::rowvec best_found_values(population_size);

  // Evaluate the initial positions.
  for (std::size_t n = 0; n < population_size; ++n) {
    const auto& parameter = positions.col(n);
    const double objective_value = problem.evaluate(parameter);
    best_found_values(n) = objective_value;

    if (objective_value <= result.objective_value) {
      result.parameter = parameter;
      result.objective_value = objective_value;
    }
  }
  ++result.iterations;
  result.duration = stopwatch.get_elapsed();

  arma::umat topology(population_size, population_size);
  bool randomize_topology = true;

  // Evaluate a single particle per loop iteration.
  while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved()) {
    if (randomize_topology) {
      topology = (arma::mat(population_size, population_size,
                            arma::fill::randu) <= neighbourhood_probability);

      // When searching for the best neighbour, we begin with the particles
      // personal best value; We don't need to visit it twice.
      topology.diag().fill(0);
      randomize_topology = false;
    }

    for (arma::uword n = 0; n < population_size; ++n) {
      // V_i^t
      const arma::vec& velocity = velocities.col(n);
      // p_i^t
      const arma::vec& best_found_parameter = best_found_parameters.col(n);
      const double best_found_value = best_found_values(n);

      // l_i^t
      arma::vec local_best_parameter = best_found_parameter;
      {
        double local_best_value = best_found_value;
        for (arma::uword i = 0; i < population_size; i++) {
          if (topology(n, i) && best_found_values(i) < local_best_value) {
            local_best_value = best_found_values(i);
            local_best_parameter = best_found_parameters.col(i);
          }
        }
      }

      const arma::vec weighted_local_attraction =
          random_uniform_in_range(0.0, maximal_local_attraction) *
          (best_found_parameter - positions.col(n));
      const arma::vec weighted_global_attraction =
          random_uniform_in_range(0.0, maximal_global_attraction) *
          (local_best_parameter - positions.col(n));
      const arma::vec acceleration =
          (weighted_local_attraction + weighted_global_attraction) / 3.0;
      const arma::vec displaced_acceleration =
          random_neighbour(acceleration, 0.0, arma::norm(acceleration));

      const double inertia = random_uniform_in_range(0, maximal_acceleration);
      velocities.col(n) = inertia * velocity + displaced_acceleration;
      positions.col(n) += velocity;
      // stay inside the bounds
      for (arma::uword k = 0; k < problem.dimension(); ++k) {
        if (positions(k, n) < problem.lower_bounds(k)) {
          positions(k, n) = problem.lower_bounds(k);
          velocities(k, n) *= -0.5;
        } else if (positions(k, n) > problem.upper_bounds(k)) {
          positions(k, n) = problem.upper_bounds(k);
          velocities(k, n) *= -0.5;
        }
      }
  
      const double objective_value = problem.evaluate(positions.col(n));
      if (objective_value < best_found_value) {
        best_found_parameters.col(n) = positions.col(n);
        best_found_values(n) = objective_value;

        if (objective_value < result.objective_value) {
          result.parameter = positions.col(n);
          result.objective_value = objective_value;
        }
      } else {
        randomize_topology = true;
      }
    }

    ++result.iterations;
    result.duration = stopwatch.get_elapsed();
  }

  return result;
}
