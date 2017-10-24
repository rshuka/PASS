#include "pass_bits/optimiser/parallel_swarm_search/parallel_swarm_search.hpp"

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

  auto start_time = std::chrono::steady_clock::now();

  // Initialisation of PSO
  pass::optimise_result result(problem.dimension());

  // Particle data, stored column-wise.
  arma::mat positions = problem.random_parameters(population_size);

  arma::mat velocities(problem.dimension(), population_size, arma::fill::randu);
  velocities *= 2 * initial_velocity;
  velocities.for_each([this](double& elem) { elem -= initial_velocity; });

  arma::mat best_found_parameters = positions;
  arma::rowvec best_found_values(population_size);

  // Evaluate the initial positions.
  for (std::size_t n = 0; n < population_size; ++n) {
    const auto& parameter = positions.col(n);
    const double objective_value = problem.evaluate(parameter);
    best_found_values(n) = objective_value;

    ++result.evaluations;

    result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start_time);

    if (objective_value <= result.objective_value) {
      result.parameter = parameter;
      result.objective_value = objective_value;

      if (result.objective_value <= acceptable_objective_value) {
        result.solved = true;
        break;
      }
    }

    if (result.evaluations >= maximal_evaluations ||
        result.iterations >= maximal_iterations ||
        result.duration >= maximal_duration) {
      break;
    }
  }
  ++result.iterations;

  // TODO: hier ändern
  while (true) {
    // MPI Call MPI_Allreduce
    for (std::size_t n = 0; n < population_size; ++n) {
      // X_i^t
      const arma::vec& position = positions.col(n);
      // V_i^t
      const arma::vec& velocity = velocities.col(n);
      // p_i^t
      const arma::vec& best_found_parameter = best_found_parameters.col(n);
      const double best_found_value = best_found_values(n);

      const arma::vec weighted_local_attraction =
          random_uniform_in_range(0.0, maximal_local_attraction) *
          (best_found_parameter - position);
      const arma::vec weighted_global_attraction =
          random_uniform_in_range(0.0, maximal_global_attraction) *
          (result.parameter - position);
      const arma::vec acceleration =
          (weighted_local_attraction + weighted_global_attraction) / 3.0;
      const arma::vec displaced_acceleration =
          random_neighbour(acceleration, 0.0, arma::norm(acceleration));

      const double inertia = random_uniform_in_range(0, maximal_acceleration);
      velocities.col(n) = inertia * velocity + displaced_acceleration;
      positions.col(n) += velocity;
      for (arma::uword k = 0; k < problem.dimension(); ++k) {
        if (position(k) < problem.lower_bounds(k)) {
          positions(k, n) = problem.lower_bounds(k);
          velocities(k, n) *= -0.5;
        } else if (position(k) > problem.upper_bounds(k)) {
          positions(k, n) = problem.upper_bounds(k);
          velocities(k, n) *= -0.5;
        }
      }

      const double objective_value = problem.evaluate(position);
      ++result.evaluations;

      result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - start_time);

      if (objective_value < best_found_value) {
        best_found_parameters.col(n) = position;
        best_found_values(n) = objective_value;
      }

      if (objective_value < result.objective_value) {
        result.parameter = position;
        result.objective_value = objective_value;
      }
    }
    if (result.objective_value <= acceptable_objective_value) {
      result.solved = true;
      break;
    }
    if (result.duration >= maximal_duration ||
        result.evaluations >= maximal_evaluations ||
        result.iterations >= maximal_iterations) {
      break;
    }
    ++result.iterations;
  }

  return result;
}