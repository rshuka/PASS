#include <pagmo/pagmo.hpp>

namespace pass
{
namespace pagmo2
{

/**
 * This class is a pagmo2 compatible UDP (user defined problem). pagmo requires
 * at least two methods of its problem type: `fitness()`, which acts the same as
 * `problem::evaluate()` in PASS, and `get_bounds()`, which returns the lower
 * and upper bounds as a tuple. This class forwards all calls to these methods
 * to `wrapped_problem`, making arbitrary PASS problems compatible to pagmo.
 *
 * This class does **not** normalise agents. Instead, the original problem
 * bounds are passed to pagmo.
 */
class problem_adapter
{
public:
  const pass::problem &wrapped_problem;

  problem_adapter(pass::problem &wrapped_problem);

  /**
   * Returns the fitness value of agent. pagmo2 supports equality and inequality
   * constraints, which would also be returned from this method call. Because
   * PASS doesn't contain any problems that use these constraints, the result
   * will always have length 1, which is the fitness value.
   */
  pagmo::vector_double fitness(const pagmo::vector_double &agent) const;

  /**
   * Returns a (lower_bounds, upper_bounds) tuple.
   */
  std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const;
}
} // namespace pagmo2
} // namespace pass
