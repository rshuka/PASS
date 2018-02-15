#pragma once

// Armadillo
#include <armadillo>

// PASS
#include "pass_bits/problem.hpp"

namespace pass {

/**
 * This problem simulates the stability of the 6PRUS robot described in
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
class parallel_kinematic_machine_6PRUS : public problem {
 public:
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
   * Caches vectors perpendicular to [redundant_joints_angles]. If you change
   * the values of [redundant_joints_angles], you *must* calculate new normals
   * and assign them to this variable!
   */
  arma::mat::fixed<3, 6> base_joints_normal;

  /**
   * Stores the lengths of the _base to middle joint_ (at index 0) and _middle
   * to platform joint_ (at index 1).
   *
   * Is initialised to (0.24, 0.56) for all three joints.
   */
  arma::mat::fixed<2, 6> link_lengths;

  /**
   * Stores the attachment points of the kinematic chains to the end effector,
   * relative to the end effector center. (`(E)x_P_i`, `(E)y_P_i`, `(E)z_P_i`)
   */
  arma::mat::fixed<3, 6> end_effector_joints_relative_position;

  /**
   * Each element stores a position in the work area that must be reached by the
   * robot as a (x, y, z, α, β, γ) vector. If any position of the trajectory
   * cannot be reached by the robot, the objective value becomes ∞.
   *
   * Is initialised to (0, 0, 0.5, 0, 0, 0).
   */
  std::vector<arma::vec::fixed<6>> end_effector_trajectory;

  /**
   * Sets the lower and upper bounds (margin for the additional actuators) to
   * (-0.6, 0.2) for all joints.
   */
  parallel_kinematic_machine_6PRUS();

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
