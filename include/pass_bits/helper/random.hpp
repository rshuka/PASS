#pragma once

#include <armadillo> // arma::vec
#include <vector>    // vector

namespace pass
{
/**
 * Returns a uniformly drawn random double number in range [min, max].
 */
double random_double_uniform_in_range(double min, double max);

/**
 * Returns a uniformly drawn random integer number in range [min, max].
 */
int random_integer_uniform_in_range(int min, int max);

/**
 * Return a randomly and uniformly drawn neighbour of *agent*, with distance
 * [*minimal_distance*, *maximal_distance*].
 */
arma::vec random_neighbour(const arma::vec &agent,
                           const double minimal_distance,
                           const double maximal_distance);

} // namespace pass
