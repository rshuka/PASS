#include "pass_bits/helper/random.hpp"
#include <cassert> // assert

std::mt19937_64 &pass::random_number_generator()
{
  static std::mt19937_64 random_number_generator;
  return random_number_generator;
}

double pass::random_uniform_in_range(double min, double max)
{
  return std::uniform_real_distribution<double>(min,
                                                max)(random_number_generator());
}

arma::vec pass::random_neighbour(const arma::vec &parameter,
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
   *  3. Translate its origin by adding *parameter*.
   */

  return parameter +
         arma::normalise(arma::vec{parameter.n_elem, arma::fill::randn}) *
             pass::random_uniform_in_range(minimal_distance, maximal_distance);
}