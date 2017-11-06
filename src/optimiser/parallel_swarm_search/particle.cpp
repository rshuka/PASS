#include "pass_bits/optimiser/parallel_swarm_search/particle.hpp"
#include "pass_bits/helper/random.hpp"

pass::particle::particle(const parallel_swarm_search& pso,
                         const pass::problem& problem)
    : pso(pso),
      problem(problem),
      position(problem.random_parameters(1)),
      velocity({problem.dimension(), arma::fill::randu}),
      best_parameter(position),
      best_value(problem.evaluate(position)) {
  velocity.for_each([&pso](double& elem) {
    elem *= 2 * pso.initial_velocity;
    elem -= pso.initial_velocity;
  });
}

bool pass::particle::update(const pass::particle& best_neighbour) {
  // p_i^t
  const arma::vec personal_weight =
      random_uniform_in_range(0.0, pso.maximal_local_attraction) *
      (best_parameter - position);

  // l_i^t
  const arma::vec neighbour_weight =
      random_uniform_in_range(0.0, pso.maximal_global_attraction) *
      (best_neighbour.best_parameter - position);

  // G_i^t
  const arma::vec gravity_center = (personal_weight + neighbour_weight) / 3.0;

  // H_i
  const arma::vec displaced_gravity_center =
      random_neighbour(gravity_center, 0.0, arma::norm(gravity_center));

  // ω
  const double inertia_weight =
      random_uniform_in_range(0, pso.maximal_acceleration);

  // V_i^{t+1}
  velocity = inertia_weight * velocity + displaced_gravity_center;
  position += velocity;

  // Check search space boundary breakouts
  for (arma::uword k = 0; k < position.n_elem; ++k) {
    if (position(k) < problem.lower_bounds(k)) {
      position(k) = problem.lower_bounds(k);
      velocity(k) = 0;
    } else if (position(k) > problem.upper_bounds(k)) {
      position(k) = problem.upper_bounds(k);
      velocity(k) = -0;
    }
  }

  const double objective_value = problem.evaluate(position);
  if (objective_value < best_value) {
    best_parameter = position;
    best_value = objective_value;
    return true;
  }
  return false;
}
