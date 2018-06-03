#include "pass_bits/helper/pagmo2/problem_adapter.hpp"
#include <stdexcept>

pass::pagmo2::problem_adapter::problem_adapter()
    : wrapped_problem(nullptr)
{
  throw std::domain_error("pass::pagmo2::problem_adapter is not really default constructible");
}

pass::pagmo2::problem_adapter::problem_adapter(const pass::problem &wrapped_problem)
    : wrapped_problem(&wrapped_problem)
{
}

pagmo::vector_double pass::pagmo2::problem_adapter::fitness(const pagmo::vector_double &agent) const
{
  return pagmo::vector_double{wrapped_problem->evaluate_normalised({agent})};
}

std::pair<pagmo::vector_double, pagmo::vector_double> pass::pagmo2::problem_adapter::get_bounds() const
{
  return {pagmo::vector_double(wrapped_problem->dimension(), 0),
          pagmo::vector_double(wrapped_problem->dimension(), 1)};
}
