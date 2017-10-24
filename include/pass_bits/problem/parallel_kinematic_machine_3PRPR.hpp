#pragma once

// Armadillo
#include <armadillo>

// PASS
#include "pass_bits/problem.hpp"

namespace pass {

/**
 * This problem simulates the stability of the 3PRPR robot described in
 * [Influence of kinematic redundancy on the singularity-free workspace of
 * parallel kinematic machines]
 * (http://link.springer.com/article/10.1007/s11465-012-0321-8).
 *
 * The underlying scenario is: The robotic arm has to move in a certain
 * trajectory, defined by any number of (x, y, γ) coordinates in the
 * 2-dimensional work area. The objective value of this problem is the
 * inaccuracy of the whole trajectory, which is defined as the maximum over the
 * inaccuracy at each intermediate point.
 *
 * Each kinematic chain is attached to an additional actuator with one degree of
 * freedom that can only be adjusted *before* the robot starts moving. Their
 * configuration influences the stability of the robot, and therefore serves as
 * the optimisation parameter.
 */
class parallel_kinematic_machine_3PRPR : public problem {
 public:
  /**
   * Stores the x/y coordinates for each base joint.
   *
   * Is initialised to (0.6, √(27)/5), (0, 0), (1.2, 0).
   */
  arma::mat::fixed<2, 3> redundant_joints_position;

  /**
   * Stores the orientation of the additional actuators. (ξ in the paper)
   *
   * Is initialised to (0, 1), (-1, 0), (-1, 0).
   */
  arma::mat::fixed<2, 3> redundant_joints_angles;

  /**
   * Stores the minimal length for each middle joint.
   *
   * Is initialised to (0.1, 0.1, 0.1).
   */
  arma::rowvec::fixed<3> middle_joints_minimal_length;

  /**
   * Stores the maximal length for each middle joint.
   *
   * Is initialised to (1.2, 1.2, 1.2).
   */
  arma::rowvec::fixed<3> middle_joints_maximal_length;

  /**
   * Stores the attachment points of the kinematic chains to the end effector,
   * relative to the end effector center.
   *
   * Is initialised to something weird, I have no idea where these values come
   * from. :(
   */
  arma::mat::fixed<2, 3> end_effector_joints_relative_position;

  /**
   * Each element stores a position in the work area that must be reached by the
   * robot as a (x, y, γ) vector, with the γ part representing the orientation
   * of the end effector. If any position of the trajectory cannot be reached by
   * the robot, the objective value becomes -∞.
   *
   * Is initialised to (0.3, 1.0, 0.0).
   */
  std::vector<arma::vec::fixed<3>> end_effector_trajectory;

  /**
   * Sets the lower and upper bounds (margin for the additional actuators) to
   * (-0.5, 0.5), (-0.2, 0.8), (-0.2, 0.8).
   */
  parallel_kinematic_machine_3PRPR();

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
