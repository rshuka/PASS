#pragma once

#include <array>
#include <mantella0>

namespace pass {
  /**
   * This optimiser works just like `mant::random_search`, but is limited to R^3
   * and uses [RFC 1149.5](https://xkcd.com/221/) to generate sample coordinates.
   *
   * Ignores `initial_parameters`.
   */
  struct rfc_1149_5_search : mant::optimiser<double, 3> {
   public:
    rfc_1149_5_search() noexcept;

//   private:
    /**
     * Holds the parameter at which problems are sampled. Is initialised with
     * the RFC 1149.5 random value in all dimensions.
     */
    static const std::array<double, 3> RANDOM_VECTOR;
  };
}
