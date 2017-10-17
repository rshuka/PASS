#include "pass_bits/problem/parallel_kinematic_machine_6PRUS.hpp"

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
// Hannover, Germany. However, we changed the coordinates of the 3rd end
// effector joint to make the configuration symmetric.
pass::parallel_kinematic_machine_6PRUS::parallel_kinematic_machine_6PRUS()
    : redundant_joints_position(
          {{-0.050659008749464, 0.050659008749464, 0.337494923062311,
            0.286835914312847, -0.286835914312847, -0.337494923062311},
           {0.360457577021932, 0.360457577021932, -0.136356800003392,
            -0.224100777018540, -0.224100777018540, -0.136356800003392},
           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}),
      redundant_joints_angles({{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                               {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                               {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}}),
      base_joints_normal(
          {{0.990268068741570, 0.990268068741570, -0.374606593415912,
            -0.615661475325659, -0.615661475325659, -0.374606593415912},
           {0.139173100960066, -0.139173100960066, -0.927183854566787,
            -0.788010753606722, 0.788010753606722, 0.927183854566787},
           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}),
      link_lengths({{0.24, 0.24, 0.24, 0.24, 0.24, 0.24},
                    {0.56, 0.56, 0.56, 0.56, 0.56, 0.56}}),
      end_effector_joints_relative_position(
          {{-0.027346281319362, 0.027346281319362, 0.072289569018135,
            0.044943287698773, -0.044943287698773, -0.072289569018135},
           {0.067684421383375, 0.067684421383375, -0.010159636370085,
            -0.057524785013291, -0.057524785013291, -0.010159636370085},
           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}),
      end_effector_trajectory({{0, 0, 0.5, 0, 0, 0}}),
      problem({-0.6, -0.6, -0.6, -0.6, -0.6, -0.6},
              {0.2, 0.2, 0.2, 0.2, 0.2, 0.2}) {}

double pass::parallel_kinematic_machine_6PRUS::evaluate(
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

    arma::mat::fixed<3, 6> middle_joints_position;
    for (arma::uword k = 0; k < base_joints_position.n_cols; ++k) {
      const std::vector<arma::vec::fixed<3>>& intersections =
          circle_sphere_intersections(
              base_joints_position.col(k), link_lengths(0, k),
              base_joints_normal.col(k), end_effector_joints_position.col(k),
              link_lengths(1, k));

      if (intersections.size() > 1) {
        middle_joints_position.col(k) = intersections.at(0);
      } else {
        return std::numeric_limits<decltype(pose_inaccuracy)>::infinity();
      }
    }

    const arma::mat::fixed<3, 6>& base_to_middle_joints_position =
        middle_joints_position - base_joints_position;
    const arma::mat::fixed<3, 6>& middle_to_end_effector_joints_position =
        end_effector_joints_position - middle_joints_position;
    const arma::mat::fixed<3, 6>& base_to_end_effector_joints_position =
        end_effector_joints_position - base_joints_position;
    const arma::mat::fixed<3, 6>& end_effector_joints_rotated_position =
        end_effector_joints_position.each_col() - end_effector_position;

    arma::mat::fixed<6, 6> forward_kinematic;
    forward_kinematic.head_rows(3) = middle_to_end_effector_joints_position;
    for (arma::uword k = 0; k < forward_kinematic.n_rows; ++k) {
      forward_kinematic.submat(3, k, 5, k) =
          arma::cross(end_effector_joints_rotated_position.col(k),
                      middle_to_end_effector_joints_position.col(k));
    }

    arma::mat::fixed<6, 12> inverse_kinematic(arma::fill::zeros);
    inverse_kinematic.diag() = base_to_end_effector_joints_position.row(0) %
                                   base_to_middle_joints_position.row(1) -
                               base_to_end_effector_joints_position.row(1) %
                                   base_to_middle_joints_position.row(0);
    for (arma::uword k = 0; k < middle_to_end_effector_joints_position.n_cols;
         ++k) {
      inverse_kinematic(k, 6 + k) =
          arma::dot(middle_to_end_effector_joints_position.col(k),
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
