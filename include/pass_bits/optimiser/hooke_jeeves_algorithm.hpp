#pragma once

#include "../optimiser.hpp"

namespace pass {
class hooke_jeeves_algorithm : public optimiser {
 public:
  /**
   * Initialized to 1.0.
   */
  double initial_stepsize;

  /**
   * Initialized to 2.0.
   */
  double stepsize_decrease;

  hooke_jeeves_algorithm() noexcept;

  virtual optimise_result optimise(const pass::problem& problem);
};
}  // namespace pass
