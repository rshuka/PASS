#pragma once

#include "pass_bits/helper/gtoc1/constants.hpp"
#include "pass_bits/problem.hpp"

namespace pass {

/**
 * `GTOC1` is an 8-dimensional optimization problem issued by the ESA:
 * http://www.esa.int/gsp/ACT/inf/projects/gtop/gtoc1.html
 * The underlying goal is to alter the trajectory of an asteroid by hitting it
 * with a spacecraft. To maximize the collision impact, the spacecraft performs
 * multiple [gravity assists](https://en.wikipedia.org/wiki/Gravity_assist).
 * This code assumes a planet sequence has been defined, and evaluates the
 * impulse gained from a whole maneuver if the spacecraft arrives at the nth
 * planet at the point in time specified by the nth problem parameter.
 */
class gtoc1 : public problem {
 public:
  /**
   * fly-by sequence of planets. This sequence is 1 element smaller than the
   * problem dimension because the last section always targets the `asteroid`.
   *
   * Is initialized to: earth, venus, earth, venus, earth, jupiter, saturn.
   */
  std::array<const celestial_body*, 7> sequence;

  /**
   * vector of flags for clockwise legs
   */
  std::array<bool, 8> rev_flag;

  /**
   * The destination of this gravity assist; the final element in the flyby
   * sequence.
   */
  asteroid destination;

  /**
   * Specific impulse of the spacecraft. Initialized to 2500.
   */
  double Isp;

  /**
   * Mass of the spacecraft. Initialized to 1500.
   */
  double mass;

  /**
   * Launch deviation of the spacecraft. Initialized to 2.5.
   */
  double DVlaunch;

  gtoc1();

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
