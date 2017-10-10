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

  pass::ackley_function problem(6);

  pass::particle_swarm_optimisation pso;
  pso.maximal_duration = std::chrono::seconds(100);
  pso.acceptable_objective_value = 1e-5;
  pso.maximal_iterations = 10;

  auto result = pso.optimise(problem);
  std::cout << "PSO found a solution of " << result.objective_value << " at "
            << result.parameter.t() << " after " << result.evaluations
            << " evaluations and"
            << " iterations " << result.iterations << std::endl;

  return 0;
}
