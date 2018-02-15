#pragma once

#include <armadillo>

namespace pass {

extern const double machine_precision;

/**
 * Rotation matrix in two dimensions
 * Further informations under:
 * (https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions)
 */
arma::mat::fixed<2, 2> rotation_matrix_2d(const double angle);

/**
 * Rotation matrix in three dimensions
 * Further informations under:
 * (https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions)
 */
arma::mat::fixed<3, 3> rotation_matrix_3d(const double roll_angle,
                                          const double pitch_angle,
                                          const double yaw_angle);

/**
 * Circle Circle Intersection
 * Further information under:
 * (http://mathworld.wolfram.com/Circle-CircleIntersection.html)
 */
std::vector<arma::vec::fixed<2>> circle_circle_intersections(
    const arma::vec::fixed<2> &first_center, const double first_radius,
    const arma::vec::fixed<2> &second_center, const double second_radius);

/**
 * Circle Sphere Intersections
 * Further informations under:
 * (https://gamedev.stackexchange.com/questions/75756/sphere-sphere-intersection-and-circle-sphere-intersection)
 */
std::vector<arma::vec::fixed<3>> circle_sphere_intersections(
    const arma::vec::fixed<3> &circle_center, const double circle_radius,
    const arma::vec::fixed<3> &circle_normal,
    const arma::vec::fixed<3> &sphere_center, const double sphere_radius);
}  // namespace pass
