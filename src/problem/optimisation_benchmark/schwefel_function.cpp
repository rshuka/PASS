#include "pass_bits/problem/optimisation_benchmark/schwefel_function.hpp"

pass::schwefel_function::schwefel_function(const arma::uword dimension)
    : problem(dimension, -500.0, 500.0, "Schwefel_Function") {}

double pass::schwefel_function::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");
  return 418.9828872724338 * dimension() -
         std::accumulate(agent.cbegin(), agent.cend(), 0.0,
                         [](const double sum, const double element) {
                           return sum + element *
                                            std::sin(std::sqrt(std::abs(element)));
                         });
}
