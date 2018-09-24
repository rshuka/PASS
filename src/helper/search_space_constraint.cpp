#include "pass_bits/helper/search_space_constraint.hpp"

#include <cmath>

arma::vec constrain_lower_bounds(const arma::vec &lower_bounds, const arma::vec &upper_bounds,
                                 arma::uword segment, arma::uword total_segments)
{
  arma::vec segment_width = (upper_bounds - lower_bounds) / arma::vec(lower_bounds.n_elem).fill(total_segments);
  arma::vec result = lower_bounds;
  result = lower_bounds + (segment - 1) * segment_width;
  return result;
}

arma::vec constrain_upper_bounds(const arma::vec &lower_bounds, const arma::vec &upper_bounds,
                                 arma::uword segment, arma::uword total_segments)
{
  arma::vec segment_width = (upper_bounds - lower_bounds) / arma::vec(lower_bounds.n_elem).fill(total_segments);
  arma::vec result = upper_bounds;
  result = lower_bounds + segment * segment_width;
  return result;
}

pass::search_space_constraint::search_space_constraint(const pass::problem &wrapped_problem,
                                                       arma::uword segment, arma::uword total_segments)
    : problem(constrain_lower_bounds(wrapped_problem.lower_bounds,
                                     wrapped_problem.upper_bounds,
                                     segment,
                                     total_segments),
              constrain_upper_bounds(wrapped_problem.lower_bounds,
                                     wrapped_problem.upper_bounds,
                                     segment,
                                     total_segments),
              wrapped_problem.name + " (segment " + std::to_string(segment) + " of " + std::to_string(total_segments) + ")"),
      wrapped_problem(wrapped_problem)
{
}

double pass::search_space_constraint::evaluate(const arma::vec &agent) const
{
  return wrapped_problem.evaluate(agent);
}
