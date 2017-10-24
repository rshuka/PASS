#include "pass_bits/problem/parallel_kinematic_machine_6PUPS.hpp"

// C++ standard library
#include <algorithm>
#include <cassert>
#include <functional>
#include <limits>
#include <string>
#include <utility>

#include "pass_bits/helper/geometry.hpp"

// The following robot configuration is based on the work of our research
// colleagues from the Institute of Mechatronic Systems, Leibniz Universität
// Hannover, Germany.
pass::parallel_kinematic_machine_6PUPS::parallel_kinematic_machine_6PUPS()
    : redundant_joints_position(
          {{-0.304055915910215, 0.304055915910215, 0.415348105304299,
            0.111292189394084, -0.111292189394084, -0.415348105304299},
           {0.304055915910215, 0.304055915910215, 0.111292189394084,
            -0.415348105304299, -0.415348105304299, 0.111292189394084},
           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}),
      redundant_joints_angles({{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                               {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                               {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}}),
      middle_joints_minimal_length({0.39, 0.39, 0.39, 0.39, 0.39, 0.39}),
      middle_joints_maximal_length({0.95, 0.95, 0.95, 0.95, 0.95, 0.95}),
      end_effector_joints_relative_position(
          {{-0.023293714059227, 0.023293714059227, 0.086933324366016,
            0.063639610306789, -0.063639610306789, -0.086933324366016},
           {0.086933324366016, 0.086933324366016, -0.023293714059227,
            -0.063639610306789, -0.063639610306789, -0.023293714059227},
           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}),
      end_effector_trajectory({{0, 0, 0.6, 0, 0, 0}}),
      problem({-0.6, -0.6, -0.6, -0.6, -0.6, -0.6},
              {0.2, 0.2, 0.2, 0.2, 0.2, 0.2}) {}

double pass::parallel_kinematic_machine_6PUPS::evaluate(
    const arma::vec& redundant_joints_actuation) const {
  assert(redundant_joints_actuation.n_elem == dimension() &&
         "`parameter` has incompatible dimension");

  double pose_inaccuracy = 0.0;
  for (const auto& end_effector_pose : end_effector_trajectory) {
    const arma::vec::fixed<3>& end_effector_position =
        end_effector_pose.head(3);
    const double end_effector_roll_angle = end_effector_pose(3);
    const double end_effector_pitch_angle = end_effector_pose(4);
    const double end_effector_yaw_angle = end_effector_pose(5);

    arma::mat::fixed<3, 6> base_joints_position = redundant_joints_position;
    for (arma::uword k = 0; k < redundant_joints_actuation.n_elem; ++k) {
      base_joints_position.col(k) +=
          redundant_joints_actuation(k) * redundant_joints_angles.col(k);
    }

    arma::mat::fixed<3, 6> end_effector_joints_position =
        rotation_matrix_3d(end_effector_roll_angle, end_effector_pitch_angle,
                           end_effector_yaw_angle) *
        end_effector_joints_relative_position;
    end_effector_joints_position.each_col() += end_effector_position;

    const arma::mat::fixed<3, 6>& base_to_end_effector_joints_position =
        end_effector_joints_position - base_joints_position;
    const arma::mat::fixed<3, 6>& end_effector_joints_rotated_position =
        end_effector_joints_position.each_col() - end_effector_position;

    const arma::rowvec::fixed<6>& middle_joints_length = arma::sqrt(
        arma::sum(arma::square(base_to_end_effector_joints_position)));
    if (arma::any(middle_joints_minimal_length > middle_joints_length ||
                  middle_joints_length > middle_joints_maximal_length)) {
      return std::numeric_limits<decltype(pose_inaccuracy)>::infinity();
    }

    arma::mat::fixed<6, 6> forward_kinematic;
    forward_kinematic.head_rows(3) = base_to_end_effector_joints_position;
    for (arma::uword k = 0; k < forward_kinematic.n_rows; ++k) {
      forward_kinematic.submat(3, k, 5, k) =
          arma::cross(end_effector_joints_rotated_position.col(k),
                      base_to_end_effector_joints_position.col(k));
    }

    arma::mat::fixed<6, 12> inverse_kinematic(arma::fill::zeros);
    inverse_kinematic.diag() = -arma::sqrt(
        arma::sum(arma::square(base_to_end_effector_joints_position)));
    for (arma::uword k = 0; k < base_to_end_effector_joints_position.n_cols;
         ++k) {
      inverse_kinematic(k, 6 + k) =
          arma::dot(base_to_end_effector_joints_position.col(k),
                    redundant_joints_angles.col(k));
    }

    arma::mat solution;
    if (!arma::solve(solution, forward_kinematic.t(), inverse_kinematic)) {
      return std::numeric_limits<decltype(pose_inaccuracy)>::infinity();
    } else {
      try {
        pose_inaccuracy = std::max(pose_inaccuracy, arma::cond(solution));
      } catch (...) {
        return std::numeric_limits<decltype(pose_inaccuracy)>::infinity();
      }
    }
  }

  return pose_inaccuracy;
}
