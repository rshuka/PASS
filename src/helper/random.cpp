#include "pass_bits/helper/random.hpp"
#include <cassert>   // assert
#include <random>    // random
#include <algorithm> // shuffle

double pass::random_double_uniform_in_range(double min, double max)
{
  return min + arma::arma_rng::randu<double>() * (max - min);
}

int pass::random_integer_uniform_in_range(int min, int max)
{
  int temp = min + arma::arma_rng::randu<double>() * (max + 1 - min);

  if (temp > max)
  {
    temp = max;
  }

  return temp;
}

arma::vec pass::random_neighbour(const arma::vec &agent,
                                 const double minimal_distance,
                                 const double maximal_distance)
{
  assert(0.0 <= minimal_distance && minimal_distance <= maximal_distance);

  /* @see J. S. Hicks and R. F. Wheeling (1959). An efficient method for
   * generating uniformly distributed points on the surface of an n-dimensional
   * sphere. Communications of the ACM, 2(4), pp. 17-19.
   *
   * In summary, it works as followed:
   *  1. Uniformly draw an *n*-dimensional direction vector (by normalising a
   *     normal distribution vector).
   *  2. Multiply it with a length, uniformly drawn
   *     from [*minimal_distance*, *maximal_distance*].
   *  3. Translate its origin by adding *agent*.
   */

  return agent +
         arma::normalise(arma::vec{agent.n_elem, arma::fill::randn}) *
             pass::random_double_uniform_in_range(minimal_distance, maximal_distance);
}
