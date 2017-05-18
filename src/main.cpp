#include <armadillo>
#include <iostream>

#include "optimiser/random_search.hpp"
#include "problem/sphere_function.hpp"

int main(int argc, char** argv) {
  pass::sphere_function problem(3);
  pass::random_search optimiser;
  optimiser.maximal_duration = std::chrono::seconds(10);

  const auto result = optimiser.optimise(problem, {});
  std::cout << "random_search found a solution of " << result.objective_value
            << " at " << result.parameter.t() << " after " << result.evaluations
            << " evaluations." << std::endl;
  return 0;
}
