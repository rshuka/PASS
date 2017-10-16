#include <armadillo>
#include <iostream>

#include "../include/pass_bits/optimiser/particle_swarm_optimisation.hpp"
#include "../include/pass_bits/optimiser/random_search.hpp"
#include "../include/pass_bits/problem/ackley_function.hpp"
#include "../include/pass_bits/problem/rastrigin_function.hpp"
#include "../include/pass_bits/problem/rosenbrock_function.hpp"
#include "../include/pass_bits/problem/sphere_function.hpp"
#include "../include/pass_bits/problem/sum_of_different_powers.hpp"

int main(int argc, char** argv) {
  arma::arma_rng::set_seed_random();

  int dimension = 6;
  pass::ackley_function problem(dimension);

  pass::particle_swarm_optimisation pso;
  pso.population_size = 40;
  pso.maximal_duration = std::chrono::seconds(10);
  pso.acceptable_objective_value = 1e-5;
  pso.maximal_evaluations = 1000000;
  int count = 0;

  // for (int n = 0; n < 51; n++) {
  auto result = pso.optimise(problem);
  std::cout << "PSO found a solution of " << result.objective_value << " at "
            << result.parameter.t() << " after " << result.evaluations
            << " evaluations." << std::endl;
  std::cout << "solved? " << result.solved << std::endl;
  if (result.objective_value <= 1e-5) {
    count++;
  }
  // }

  std::cout << "Count " << count << std::endl;
  return 0;
}
