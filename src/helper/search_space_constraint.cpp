#include "pass_bits/helper/search_space_constraint.hpp"

#include <cmath>

arma::vec constrain_lower_bounds(const arma::vec &lower_bounds, const arma::vec &upper_bounds,
                                 arma::uword segment, arma::uword total_segments)
{
  double segment_width = std::abs(upper_bounds[0] - lower_bounds[0]) / total_segments;
  arma::vec result = lower_bounds;
  result[0] = lower_bounds[0] + (segment - 1) * segment_width;
  return result;
}

arma::vec constrain_upper_bounds(const arma::vec &lower_bounds, const arma::vec &upper_bounds,
                                 arma::uword segment, arma::uword total_segments)
{
  double segment_width = std::abs(upper_bounds[0] - lower_bounds[0]) / total_segments;
  arma::vec result = upper_bounds;
  result[0] = lower_bounds[0] + segment * segment_width;
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
