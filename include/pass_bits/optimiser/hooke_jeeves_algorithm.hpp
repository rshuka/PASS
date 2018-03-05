#pragma once

#include "pass_bits/optimiser.hpp"

namespace pass
{
/**
 * Implements the Pattern Search (Hooke Jeeves) algorithm. 
 * (https://en.wikipedia.org/wiki/Pattern_search_(optimization))
 * The initial stepsize depends on the problem boundarys and 
 * the stepsize decreases by half if no better solution 
 * is found.
 */  
class hooke_jeeves_algorithm : public optimiser
{
public:
  virtual optimise_result optimise(const pass::problem &problem);
};
} // namespace pass
