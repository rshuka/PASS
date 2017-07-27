#pragma once

#include "../problem.hpp"

namespace pass {

/**
 * This problem simulates the stability of the 3PRRR robot described in
 * [Influence of kinematic redundancy on the singularity-free workspace of
 * parallel kinematic machines]
 * (http://link.springer.com/article/10.1007/s11465-012-0321-8).
 */
class parallel_kinematic_machine_3PRRR : public problem {
 public:
  arma::mat::fixed<2, 3> redundant_joints_position;
  arma::mat::fixed<2, 3> redundant_joints_angles;
  arma::mat::fixed<2, 3> link_lengths;
  arma::mat::fixed<2, 3> end_effector_joints_relative_position;
  std::vector<arma::vec::fixed<3>> end_effector_trajectory;

  parallel_kinematic_machine_3PRRR();

  virtual double evaluate(const arma::vec& parameter) const override;
};
}  // namespace pass
