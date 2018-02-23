#pragma once

#include "pass_bits/optimiser.hpp"

namespace pass
{
/**
 * Implements the Random Search Algorithm
 * (https://en.wikipedia.org/wiki/Random_search)
 */
class random_search : public optimiser
{
public:
  virtual optimise_result optimise(const pass::problem &problem);
};
} // namespace pass
