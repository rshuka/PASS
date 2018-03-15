#include "pass_bits/problem/rosenbrock_function.hpp"

pass::rosenbrock_function::rosenbrock_function(const arma::uword dimension)
    : problem(dimension, -2.048, 2.048) {}

double pass::rosenbrock_function::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");
  return std::inner_product(
      agent.cbegin(), std::prev(agent.cend(), 1),
      std::next(agent.cbegin(), 1), 0.0, std::plus<double>(),
      [](const double element, const double other_element) {
        return 100.0 * std::pow(other_element - std::pow(element, 2.0), 2.0) +
               std::pow(element - 1.0, 2.0);
      });
}
