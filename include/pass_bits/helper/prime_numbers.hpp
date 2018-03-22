#pragma once

#include <array>
#include <armadillo>

namespace pass
{
/**
 * Primes for use in `pass::problem::hammersley_agents`. Taken from
 * https://cis.temple.edu/~beigel/cis573/alizzy/prime-list.html
 */
extern const std::array<arma::uword, 2000> prime_numbers;

} // namespace pass
