#pragma once

#include <random>    // std::mt19937_64
#include <armadillo> // arma::sed

namespace pass
{
class seed
{
public:
  // The random number generator should act as a singleton, no need for default constructors.
  seed() = delete;
  seed(const seed &) = delete;
  seed &operator=(const seed &) = delete;

  // The std::mersenne_twister_engine requires an l-value (a reference in this case) instead of an r-value.
  // Otherwise, calls like std::uniform_real_distribution<double>(0, 1)(seed::get_generator()) won't compile.
  static std::mt19937_64 &get_generator();

  static void set_seed(
      const arma::arma_rng::seed_type seed);

  static void set_random_seed();

  static arma::arma_rng::seed_type get_seed();

protected:
  static arma::arma_rng::seed_type seed_;
  static std::mt19937_64 generator_;
};
} // namespace pass
