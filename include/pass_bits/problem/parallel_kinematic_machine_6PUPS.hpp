#pragma once

// Armadillo
#include <armadillo>

// PASS
#include "pass_bits/problem.hpp"

namespace pass {

/**
 * This problem simulates the stability of the 6PUPS robot described in
 * [Influence of kinematic redundancy on the singularity-free workspace of
 * parallel kinematic machines]
 * (http://link.springer.com/article/10.1007/s11465-012-0321-8).
 *
 * The underlying scenario is: The robotic arm has to move in a certain
 * trajectory, defined by any number of (x, y, z, α, β, γ) coordinates in the
 * 3-dimensional work area. The objective value of this problem is the
 * inaccuracy of the whole trajectory, which is defined as the maximum over the
 * inaccuracy at each intermediate point.
 *
 * Each kinematic chain is attached to an additional actuator with one degree of
 * freedom that can only be adjusted *before* the robot starts moving. Their
 * configuration influences the stability of the robot, and therefore serves as
 * the optimisation parameter.
 */
class parallel_kinematic_machine_6PUPS : public problem {
 public:
  /**
   * Initialises this object according to table 5 in the paper.
   */
  parallel_kinematic_machine_6PUPS();

  /**
   * Stores the x/y/z coordinates for each base joint. (`(0)x_G_i`, `(0)y_G_i`,
   * `(0)z_G_i`)
   */
  arma::mat::fixed<3, 6> redundant_joints_position;

  /**
   * Stores the orientation of the additional actuators. (ξ and χ; the third
   * element (orientation in z-axis) is unnamed in the paper)
   */
  arma::mat::fixed<3, 6> redundant_joints_angles;

  /**
   * Stores the minimum lengths of each actuator. (ρ_min)
   *
   * Is initialised to 0.39 for all joints.
   */
  arma::rowvec::fixed<6> middle_joints_minimal_length;

  /**
   * Stores the maximum lengths of each actuator. (ρ_max)
   *
   * Is initialised to 0.95 for all joints.
   */
  arma::rowvec::fixed<6> middle_joints_maximal_length;

  /**
   * Stores the attachment points of the kinematic chains to the end effector,
   * relative to the end effector center. (`(E)x_P_i`, `(E)y_P_i`, `(E)z_P_i`)
   */
  arma::mat::fixed<3, 6> end_effector_joints_relative_position;

  /**
   * Each element stores a position in the work area that must be reached by the
   * robot as a (x, y, z, α, β, γ) vector. If any position of the trajectory
   * cannot be reached by the robot, the objective value becomes -∞.
   *
   * Is initialised to (0, 0, 0.6, 0, 0, 0).
   */
  std::vector<arma::vec::fixed<6>> end_effector_trajectory;

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
