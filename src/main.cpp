#include <iostream>
#include <mantella0>
#include "optimiser/rfc_1149_5.hpp"

int main() {
  mant::sphere_function<double, 3> problem;
  pass::rfc_1149_5_search solver;
  std::array<double, 3> initial_parameter = {0, 0, 0};
  const auto result =
      solver.optimisation_function(problem, {initial_parameter});

  std::cout << "Found an objective value of " << result.objective_value
            << " at " << result.parameter[0] << ", " << result.parameter[1]
            << ", " << result.parameter[2] << "." << std::endl;
  return 0;
}
