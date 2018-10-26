#include "pass_bits/helper/random.hpp"
#include <cassert>   // assert
#include <random>    // random
#include <algorithm> // shuffle

double pass::random_double_uniform_in_range(double min, double max)
{
  assert(min < max && "'min' should be less than 'max'");

  return min + arma::arma_rng::randu<double>() * (max - min);
}

int pass::random_integer_uniform_in_range(int min, int max)
{
  assert(min < max && "'min' should be less than 'max'");

  int temp = min + arma::arma_rng::randu<double>() * (max + 1 - min);

  if (temp > max)
  {
    temp = max;
  }

  assert(temp <= max && temp >= min && "Value not in range");

  return temp;
}

arma::rowvec pass::integers_uniform_in_range(const int min, const int max, const int count)
{
  assert(min < max && "'min' should be less than 'max'");
  assert(count > 1 && "'count' should be more than 1");

  arma::rowvec numbers(count);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(min, max);

  for (int n = 0; n < count; ++n)
  {
    numbers(n) = dis(gen);
  }

  return numbers;
}

arma::vec pass::random_neighbour(const arma::vec &agent,
                                 const double minimal_distance,
                                 const double maximal_distance)
{
  assert(0.0 <= minimal_distance && minimal_distance <= maximal_distance && "Neighbour not in the range");

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
