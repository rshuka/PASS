#include "../../include/pass_bits/optimiser/particle_swarm_optimisation.hpp"

// pass::random_number_generator(), pass::random_neighbour()
#include "../../include/pass_bits/helper/random.hpp"

// std::accumulate
#include <algorithm>

// std::pow
#include <cmath>

// assert
#include <cassert>

pass::particle_swarm_optimisation::particle_swarm_optimisation() noexcept
    : optimiser(),
      initial_velocity(0.5),
      maximal_acceleration(1.0 / (2.0 * std::log(2.0))),
      maximal_local_attraction(0.5 + std::log(2.0)),
      maximal_global_attraction(maximal_local_attraction) {}

pass::optimise_result pass::particle_swarm_optimisation::optimise(
    const pass::problem& problem, const arma::mat& initial_parameters) {
  const arma::uword dimension = problem.dimension();
  const arma::uword particle_count = initial_parameters.n_cols;

  assert(initial_velocity >= 0.0);
  assert(maximal_acceleration >= 0.0);
  assert(maximal_local_attraction >= 0.0);
  assert(maximal_global_attraction >= 0.0);
  assert(initial_parameters.n_rows == dimension &&
         initial_parameters.n_cols > 0 &&
         "`initial_parameters` must have the same dimension as `problem` and "
         "not be empty");

  auto start_time = std::chrono::steady_clock::now();
  pass::optimise_result result(dimension);

  // Best found parameters and objective values for each particle. The data for
  // each particle is stored in a column.
  arma::mat best_found_parameters = initial_parameters;
  arma::rowvec best_found_values(particle_count);

  // The current positions and velocities of each particle.
  arma::mat positions = best_found_parameters;
  arma::mat velocities{arma::size(best_found_parameters), arma::fill::randu};
  velocities *= 2 * initial_velocity;
  velocities.for_each([this](auto& elem) { elem -= initial_velocity; });

  // Evaluate the initial parameters.
  for (arma::uword n = 0; n < particle_count; ++n) {
    const auto& parameter = best_found_parameters.col(n);
    const double objective_value = best_found_values(n) =
        problem.evaluate(parameter);

    ++result.evaluations;
    result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start_time);

    if (objective_value <= result.objective_value) {
      result.parameter = parameter;
      result.objective_value = objective_value;

      if (result.objective_value <= acceptable_objective_value) {
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

  // Evaluate a single particle per iteration.
  while (result.duration < maximal_duration &&
         result.evaluations < maximal_evaluations &&
         result.objective_value > acceptable_objective_value) {
    const auto n = result.evaluations % particle_count;
    const auto& position = positions.col(n);
    const auto& velocity = velocities.col(n);
    const auto& best_found_parameter = best_found_parameters.col(n);
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
    for (std::size_t k = 0; k < dimension; ++k) {
      if (position(k) < problem.lower_bounds(k)) {
        positions(n, k) = problem.lower_bounds(k);
        velocities(n, k) *= -0.5;
      } else if (position(k) > problem.upper_bounds(k)) {
        positions(n, k) = problem.upper_bounds(k);
        velocities(n, k) *= -0.5;
      }
    }

    const double objective_value = problem.evaluate(position);
    ++result.evaluations;
    if (n == particle_count - 1) {
      ++result.iterations;
    }
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

  return result;
}
