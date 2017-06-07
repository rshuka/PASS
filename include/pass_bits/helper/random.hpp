// std::mt19937_64
#include <random>

// arma::vec
#include <armadillo>

namespace pass {
/**
 * Returns a reference to a single random number generator.
 *
 * TODO: Right now, some functions use `arma::arma_rng` instead for random
 * numbers.
 *
 * TODO: Check everything for thread safety when we start parallelising.
 */
inline std::mt19937_64& random_number_generator();

/**
 * Return a randomly and uniformly drawn neighbour of *parameter*, with distance
 * [*minimal_distance*, *maximal_distance*].
 */
arma::vec random_neighbour(const arma::vec& parameter,
                           const double minimal_distance,
                           const double maximal_distance);
}  // namespace pass
