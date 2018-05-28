#include "pass_bits/helper/pagmo2/problem_adapter.hpp"
#include <armadillo>

pass::pagmo2::problem_adapter::problem_adapter(const pass::problem &wrapped_problem)
    : wrapped_problem(wrapped_problem);

pagmo::vector_double pass::pagmo2::problem_adapter::fitness(const pagmo::vector_double &agent) const
{
  return wrapped_problem.evaluate({agent});
}

std::pair<pagmo::vector_double, pagmo::vector_double> pass::pagmo2::problem_adapter::get_bounds() const
{
  return {
    arma::conv_to<pagmo::vector_double>(lower_bounds),
    arma::conv_to<pagmo::vector_double>(upper_bounds)
  };
}
