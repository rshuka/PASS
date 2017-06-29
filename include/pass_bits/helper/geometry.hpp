#include <armadillo>

namespace pass {
extern const double machine_precision;

arma::mat::fixed<2, 2> rotation_matrix_2d(const double angle);

arma::mat::fixed<3, 3> rotation_matrix_3d(const double roll_angle,
                                          const double pitch_angle,
                                          const double yaw_angle);

std::vector<arma::vec::fixed<2>> circle_circle_intersections(
    const arma::vec::fixed<2>& first_center, const double first_radius,
    const arma::vec::fixed<2>& second_center, const double second_radius);

std::vector<arma::vec::fixed<3>> circle_sphere_intersections(
    const arma::vec::fixed<3>& circle_center, const double circle_radius,
    const arma::vec::fixed<3>& circle_normal,
    const arma::vec::fixed<3>& sphere_center, const double sphere_radius);
}  // namespace pass
