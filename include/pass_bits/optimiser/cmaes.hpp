#pragma once

#include "pass_bits/optimiser.hpp"

namespace pass
{
/**
 *
 */
class cmaes : public optimiser
{
public:
  /**
   * Must be at least 5. Is initialized to 0.
   */
  arma::uword population_size;

  /**
   * Initialises the optimiser with its name.
   */
  cmaes() noexcept;

  virtual optimise_result optimise(const pass::problem &problem);
};
} // namespace pass
