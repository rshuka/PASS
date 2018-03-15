#pragma once

#include <armadillo> // arma::vec

namespace pass
{
/**
 * Returns a uniformly drawn random number in range [min, max].
 */
double random_uniform_in_range(double min, double max);

/**
 * Return a randomly and uniformly drawn neighbour of *agent*, with distance
 * [*minimal_distance*, *maximal_distance*].
 */
arma::vec random_neighbour(const arma::vec &agent,
                           const double minimal_distance,
                           const double maximal_distance);

} // namespace pass
