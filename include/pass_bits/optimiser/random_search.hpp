#pragma once

#include "../optimiser.hpp"

namespace pass {
class random_search : public optimiser {
 public:
  virtual optimise_result optimise(const pass::problem& problem,
                                   const arma::mat& initial_parameters);
};
}  // namespace pass
