#include <armadillo>
#include <iostream>

#include "../include/pass_bits/optimiser/particle_swarm_optimisation.hpp"
#include "../include/pass_bits/optimiser/random_search.hpp"
#include "../include/pass_bits/problem/sphere_function.hpp"

int main(int argc, char** argv) {
  arma::arma_rng::set_seed_random();

  pass::sphere_function problem(3);

  pass::random_search random_search;
  random_search.maximal_duration = std::chrono::seconds(10);
  auto result = random_search.optimise(problem, {});
  std::cout << "random_search found a solution of " << result.objective_value
            << " at " << result.parameter.t() << " after " << result.evaluations
            << " evaluations." << std::endl;

  pass::particle_swarm_optimisation pso;
  pso.maximal_duration = std::chrono::seconds(10);
  result = pso.optimise(problem,
                        problem.random_parameters(problem.dimension() * 20));
  std::cout << "PSO found a solution of " << result.objective_value << " at "
            << result.parameter.t() << " after " << result.evaluations
            << " evaluations." << std::endl;

  return 0;
}
