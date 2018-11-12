#include "pass_bits/helper/seed.hpp"

namespace pass
{
decltype(seed::seed_) seed::seed_ = 12345;
decltype(seed::generator_) seed::generator_;

std::mt19937_64 &seed::get_generator()
{
  return generator_;
}

void seed::set_seed(const arma::arma_rng::seed_type seed)
{
  seed_ = seed;

  generator_.seed(seed_);
  arma::arma_rng::set_seed(seed_);
}

void seed::set_random_seed()
{

  arma::arma_rng::set_seed_random();

#if defined(SUPPORT_MPI)
  set_seed(arma::randi<arma::Col<arma::arma_rng::seed_type>>(static_cast<arma::uword>(pass::number_of_nodes()))(static_cast<arma::uword>(pass::node_rank())));
#else
  set_seed(arma::randi<arma::Col<arma::arma_rng::seed_type>>(1)(0));
#endif
}

arma::arma_rng::seed_type seed::get_seed()
{
  return seed_;
}
} // namespace pass
