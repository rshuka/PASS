#include "pass_bits/optimiser/parallel_swarm_search/particle.hpp"
#include "pass_bits/helper/random.hpp"

pass::particle::particle(const parallel_swarm_search &pss,
                         const pass::problem &problem)
    : pss(pss),
      problem(problem),
      position(problem.random_agents(1)),
      velocity(arma::vec(problem.dimension())),
      personal_best_position(position),
      personal_best_fitness_value(problem.evaluate(position))
{
  for (arma::uword i = 0; i < problem.dimension(); i++)
  {
    velocity[i] = random_uniform_in_range(
        problem.lower_bounds[i] - position[i],
        problem.upper_bounds[i] - position[i]);
  }
}

bool pass::particle::update(const pass::particle &local_best_agent)
{
  //p_i
  const arma::vec weighted_personal_attraction = position +
                                                 random_uniform_in_range(0.0, pss.cognitive_acceleration) *
                                                     (personal_best_position - position);

  // l_i
  const arma::vec weighted_local_attraction = position +
                                              random_uniform_in_range(0.0, pss.social_acceleration) *
                                                  (local_best_agent.personal_best_position - position);

  // G
  arma::vec attraction_center;

  // If the best informant is the particle itself, define the gravity center G as the middle of x-p'
  if (&local_best_agent == this)
  {
    attraction_center = 0.5 * (position + weighted_personal_attraction);
  }
  else
  {
    attraction_center = (position + weighted_personal_attraction + weighted_local_attraction) / 3.0;
  }

  velocity = pss.inertia * velocity +
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

  const double fitness_value = problem.evaluate(position);
  if (fitness_value < personal_best_fitness_value)
  {
    personal_best_position = position;
    personal_best_fitness_value = fitness_value;
    return true;
  }
  return false;
}
