#include "pass_bits/problem.hpp"
#include "pass_bits/helper/random.hpp"
#include "pass_bits/helper/prime_numbers.hpp"

arma::vec pass::problem::bounds_range() const noexcept
{
  return upper_bounds - lower_bounds;
}

arma::uword pass::problem::dimension() const noexcept
{
  return lower_bounds.n_elem;
}

pass::problem::problem(const arma::uword dimension, const double lower_bound,
                       const double upper_bound, const std::string name)
    : lower_bounds(arma::vec(dimension).fill(lower_bound)),
      upper_bounds(arma::vec(dimension).fill(upper_bound)),
      name(name)
{
  assert(lower_bound < upper_bound &&
         "`problem.lower_bounds` must be less than `problem.upper_bounds`");
}

pass::problem::problem(const arma::vec &lower_bounds,
                       const arma::vec &upper_bounds, const std::string name)
    : lower_bounds(lower_bounds),
      upper_bounds(upper_bounds),
      name(name)
{
  assert(upper_bounds.n_elem == dimension() &&
         "`problems.lower_bounds` and `problem.upper_bounds` must have the "
         "same dimension");
  assert(arma::all(lower_bounds < upper_bounds) &&
         "`problem.lower_bounds` must be less than `problem.upper_bounds` in "
         "each dimension");
}

double pass::problem::evaluate_normalised(const arma::vec &agent) const
{
  return evaluate(agent % bounds_range() + lower_bounds);
}

arma::mat pass::problem::normalised_random_agents(const arma::uword count) const
{
  assert(count >= 1 && "Can't generate 0 agents");
  return arma::mat{dimension(), count, arma::fill::randu};
}

arma::mat pass::problem::normalised_hammersley_agents(const arma::uword count) const
{
  assert(count >= 1 && "Can't generate 0 agents");
  assert(dimension() + 1 <= pass::prime_numbers.size() &&
         "Can't generate hammersley points because `dimension` "
         "exceeds the number of available primes");

  // Calculate hammersley points.
  arma::mat agents(dimension(), count);

  // The first dimension is in the pure Hammersley particular
  // A random number define which dimension is "The First"
  arma::uword random_dimension = pass::random_integer_uniform_in_range(0, agents.n_rows - 1);

  for (arma::uword k = 0; k < count; k++)
  {
    agents(random_dimension, k) = static_cast<double>(k) / static_cast<double>(count);

    for (arma::uword d = 0; d < agents.n_rows; d++)
    {
      if (d == random_dimension)
      {
        continue;
      }
      arma::uword p = pass::prime_numbers[pass::random_integer_uniform_in_range(0, agents.n_rows - 1)];
      double phi = 0;

      // Calculate Φ_p(k) - Pseudocode on page 3
      // Compare http://www.cse.cuhk.edu.hk/~ttwong/papers/udpoint/udpoint.pdf,
      {
        arma::uword p_temp = p;
        arma::uword k_temp = k;
        while (k_temp > 0)
        {
          phi += static_cast<double>(k_temp % p) / p_temp;
          k_temp /= p;
          p_temp *= p;
        }
      }

      agents(d, k) = phi;
    }
  }

  return agents;
}

arma::mat pass::problem::initialise_normalised_agents(const arma::uword count) const
{
  if (arma::randu() <= 0.5)
  {
    return normalised_hammersley_agents(count);
  }
  else
  {
    return normalised_random_agents(count);
  }
}
