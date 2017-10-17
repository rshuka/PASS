#include "pass_bits/problem/parallel_kinematic_machine_3PRRR.hpp"

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
pass::parallel_kinematic_machine_3PRRR::parallel_kinematic_machine_3PRRR()
    : redundant_joints_position(
          {{0.6, 0.0, 1.2}, {1.039230484541327, 0.0, 0.0}}),
      redundant_joints_angles({{0.0, -1.0, -1.0}, {1.0, 0.0, 0.0}}),
      link_lengths({{0.6, 0.6, 0.6}, {0.6, 0.6, 0.6}}),
      end_effector_joints_relative_position(
          {{-0.000066580445834, -0.092751709777083, 0.092818290222917},
           {0.106954081945581, -0.053477040972790, -0.053477040972790}}),
      end_effector_trajectory({{0.3, 1.0, 0.0}}),
      problem({-0.5, -0.2, -0.2}, {0.5, 0.8, 0.8}) {}

double pass::parallel_kinematic_machine_3PRRR::evaluate(
    const arma::vec& redundant_joints_actuation) const {
  assert(redundant_joints_actuation.n_elem == dimension() &&
         "`parameter` has incompatible dimension");

  double pose_inaccuracy = 0.0;
  for (const auto& end_effector_pose : end_effector_trajectory) {
    const arma::vec::fixed<2>& end_effector_position =
        end_effector_pose.head(2);
    const double end_effector_angle = end_effector_pose(2);

    arma::mat::fixed<2, 3> base_joints_position = redundant_joints_position;
    for (arma::uword k = 0; k < redundant_joints_actuation.n_elem; ++k) {
      base_joints_position.col(k) +=
          redundant_joints_actuation(k) * redundant_joints_angles.col(k);
    }

    arma::mat::fixed<2, 3> end_effector_joints_position =
        rotation_matrix_2d(end_effector_angle) *
        end_effector_joints_relative_position;
    end_effector_joints_position.each_col() += end_effector_position;

    arma::mat::fixed<2, 3> middle_joints_position;
    for (arma::uword k = 0; k < base_joints_position.n_cols; ++k) {
      const std::vector<arma::vec::fixed<2>>& intersections =
          circle_circle_intersections(
              base_joints_position.col(k), link_lengths(0, k),
              end_effector_joints_position.col(k), link_lengths(1, k));

      if (intersections.size() > 1) {
        middle_joints_position.col(k) = intersections.at(0);
      } else {
        return std::numeric_limits<decltype(pose_inaccuracy)>::infinity();
      }
    }

    const arma::mat::fixed<2, 3>& base_to_middle_joints_position =
        middle_joints_position - base_joints_position;
    const arma::mat::fixed<2, 3>& middle_to_end_effector_joints_position =
        end_effector_joints_position - middle_joints_position;
    const arma::mat::fixed<2, 3>& end_effector_joints_rotated_position =
        end_effector_joints_position.each_col() - end_effector_position;

    arma::mat::fixed<3, 3> forward_kinematic;
    forward_kinematic.head_rows(2) = middle_to_end_effector_joints_position;
    forward_kinematic.row(2) =
        -forward_kinematic.row(0) %
            end_effector_joints_rotated_position.row(1) +
        forward_kinematic.row(1) % end_effector_joints_rotated_position.row(0);

    arma::mat::fixed<3, 6> inverse_kinematic(arma::fill::zeros);
    inverse_kinematic.diag() =
        forward_kinematic.row(0) % base_to_middle_joints_position.row(1) -
        forward_kinematic.row(1) % base_to_middle_joints_position.row(0);
    for (arma::uword k = 0; k < redundant_joints_actuation.n_elem; ++k) {
      inverse_kinematic(k, 3 + k) =
          -arma::dot(middle_to_end_effector_joints_position.col(k),
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
