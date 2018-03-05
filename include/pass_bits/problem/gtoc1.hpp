#pragma once

#include "pass_bits/helper/gtoc1/constants.hpp"
#include "pass_bits/problem.hpp"

namespace pass
{
/**
 * `GTOC1` is an 8-dimensional optimization problem issued by the ESA:
 * http://www.esa.int/gsp/ACT/inf/projects/gtop/gtoc1.html
 * The underlying goal is to alter the trajectory of an asteroid by hitting it
 * with a spacecraft. To maximize the collision impact, the spacecraft performs
 * multiple gravity assists (https://en.wikipedia.org/wiki/Gravity_assist).
 * This code assumes a planet sequence has been defined, and evaluates the
 * impulse gained from a whole maneuver if the spacecraft arrives at the nth
 * planet at the point in time specified by the nth problem parameter.
 */
class gtoc1 : public problem
{
public:
  /**
   * Fly-by sequence of planets. This sequence is 1 element smaller than the
   * problem dimension because the last section always targets the `asteroid`.
   *
   * Is initialized to: earth, venus, earth, venus, earth, jupiter, saturn.
   */
  std::array<const celestial_body *, 7> sequence;

  /**
   * Vector of flags for clockwise legs
   * Tell if the planet is passed in a clockwise way or not
   * 1 - clockwise
   * 0 - counter clockwise
   */
  std::array<bool, 8> rev_flag;

  /**
   * The destination of this gravity assist; the final element in the flyby
   * sequence.
   * (First Vector) Data for the Asteroid 2001 TW229
   * (Second Double Value) Asteroid Epoche initialized with 53600.0
   * (Third Double Value) mu (standard gravitational parameter) initialized with
   * 0.0 Further Info:
   * http://www.esa.int/gsp/ACT/doc/MAD/ACT-MEM-MAD-GTOC1-The%20Problem_V4.pdf
   */
  asteroid destination;

  /**
   * Specific impulse of the spacecraft. Initialized to 2500 Sek.
   */
  double Isp;

  /**
   * Mass of the spacecraft. Initialized to 1500 Kg.
   */
  double mass;

  /**
   * Launch deviation of the spacecraft. Initialized to 2.5 km/sec.
   */
  double DVlaunch;

  /**
   * Sets the lower and upper bounds for the times to
   * (3000, 10000) (14, 2000) (14, 2000) (14, 2000)
   * (14, 2000) (100, 9000) (366, 9000) (300, 9000)
   *
   * t0: start time (in MJD2000)
   * t1-7: transfer time from planet to another planet (in days)
   */
  gtoc1();

  virtual double evaluate(const arma::vec &agent) const override;
};
} // namespace pass
