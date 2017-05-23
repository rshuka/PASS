#include <armadillo>
#include <iostream>

#include "../include/pass_bits/optimiser/random_search.hpp"
#include "../include/pass_bits/problem/sphere_function.hpp"

int main(int argc, char** argv) {
  arma::arma_rng::set_seed_random();

  pass::sphere_function problem(3);
  pass::random_search optimiser;

  const auto result = optimiser.optimise(problem, {});
  std::cout << "random_search found a solution of " << result.objective_value
            << " at " << result.parameter.t() << " after " << result.evaluations
            << " evaluations." << std::endl;
  return 0;
}
