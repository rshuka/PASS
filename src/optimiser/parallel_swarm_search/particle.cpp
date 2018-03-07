#include "pass_bits/optimiser/parallel_swarm_search/particle.hpp"
#include "pass_bits/helper/random.hpp"

pass::particle::particle(const parallel_swarm_search &pso,
                         const pass::problem &problem)
    : pso(pso),
      problem(problem),
      position(problem.random_agents(1)),
      velocity(arma::vec(problem.dimension())),
      best_agent(position),
      best_value(problem.evaluate(position))
{
  for (arma::uword i = 0; i < problem.dimension(); i++) {
    velocity[i] = random_uniform_in_range(
        problem.lower_bounds[i] - position[i],
        problem.upper_bounds[i] - position[i]);
  }
}

bool pass::particle::update(const pass::particle &best_neighbour)
{
  //p_i
  const arma::vec weighted_personal_attraction = position +
      random_uniform_in_range(0.0, pso.cognitive_acceleration) *
        (best_agent - position);

  // l_i
  const arma::vec weighted_local_attraction = position +
      random_uniform_in_range(0.0, pso.social_acceleration) *
        (best_neighbour.best_agent - position);

  // G
  arma::vec attraction_center;

  // If the best informant is the particle itself, define the gravity center G as the middle of x-p'
  if (&best_neighbour == this)
  {
    attraction_center = 0.5 * (position + weighted_personal_attraction);
  }
  else
  {
    attraction_center = (position + weighted_personal_attraction + weighted_local_attraction) / 3.0;
  }

  velocity = pso.inertia * velocity +
             random_neighbour(attraction_center, 0.0, arma::norm(attraction_center - position)) -
             position;

  // move by applying this new velocity to the current position
  position += velocity;

  // Check search space boundary breakouts
  for (arma::uword k = 0; k < position.n_elem; ++k)
  {
    if (position(k) < problem.lower_bounds(k))
    {
      position(k) = problem.lower_bounds(k);
      velocity(k) *= -0.5;
    }
    else if (position(k) > problem.upper_bounds(k))
    {
      position(k) = problem.upper_bounds(k);
      velocity(k) *= -0.5;
    }
  }

  const double objective_value = problem.evaluate(position);
  if (objective_value < best_value)
  {
    best_agent = position;
    best_value = objective_value;
    return true;
  }
  return false;
}
