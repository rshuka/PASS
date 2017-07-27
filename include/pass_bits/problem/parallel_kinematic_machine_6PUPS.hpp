#pragma once

#include "../problem.hpp"

namespace pass {

/**
 * This problem simulates the stability of the 6PUPS robot described in
 * [Influence of kinematic redundancy on the singularity-free workspace of
 * parallel kinematic machines]
 * (http://link.springer.com/article/10.1007/s11465-012-0321-8).
 */
class parallel_kinematic_machine_6PUPS : public problem {
 public:
  arma::mat::fixed<3, 6> redundant_joints_position;
  arma::mat::fixed<3, 6> redundant_joints_angles;
  arma::rowvec::fixed<6> middle_joints_minimal_length;
  arma::rowvec::fixed<6> middle_joints_maximal_length;
  arma::mat::fixed<3, 6> end_effector_joints_relative_position;
  std::vector<arma::vec::fixed<6>> end_effector_trajectory;

  parallel_kinematic_machine_6PUPS();

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
