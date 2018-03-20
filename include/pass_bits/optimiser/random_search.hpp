#pragma once

#include "pass_bits/optimiser.hpp"

namespace pass
{
/**
 * Implements the Random Search algorithm
 * (https://en.wikipedia.org/wiki/Random_search)
 */
class random_search : public optimiser
{
public:
  /**
   * Initialises the optimiser with its name
   */
  random_search() noexcept;

  virtual optimise_result optimise(const pass::problem &problem);
};
} // namespace pass
