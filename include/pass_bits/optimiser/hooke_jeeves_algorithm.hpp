#pragma once

#include "pass_bits/optimiser.hpp"

namespace pass {
class hooke_jeeves_algorithm : public optimiser {
 public:
  /**
   * Initialized to 0.0. Should be set to `problem.bounds_range() / 2` by the
   * caller.
   */
  double initial_stepsize;

  /**
   * In iterations without improvement in any direction, divide the stepsize by
   * this value.
   *
   * Initialized to 2.0.
   */
  double stepsize_decrease;

  hooke_jeeves_algorithm() noexcept;

  virtual optimise_result optimise(const pass::problem& problem);
};
}  // namespace pass
