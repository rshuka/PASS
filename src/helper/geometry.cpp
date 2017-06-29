#include "../../include/pass_bits/helper/geometry.hpp"

// C++ standard library
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace pass {
const double machine_precision = 1e-12;

arma::mat::fixed<2, 2> rotation_matrix_2d(const double angle) {
  return {{std::cos(angle), -std::sin(angle)},
          {std::sin(angle), std::cos(angle)}};
}

arma::mat::fixed<3, 3> rotation_matrix_3d(const double roll_angle,
                                          const double pitch_angle,
                                          const double yaw_angle) {
  assert(std::isfinite(roll_angle) &&
         "rotation_matrix_3d: The roll angle must be finite.");
  assert(std::isfinite(pitch_angle) &&
         "rotation_matrix_3d: The pitch angle must be finite.");
  assert(std::isfinite(yaw_angle) &&
         "rotation_matrix_3d: The yaw angle must be finite.");

  // In case the Tait-Bryan angles losses a rank, i.e. a gimbal lock would
  // occur.
  if (std::abs(std::fmod(pitch_angle, arma::datum::pi / 2.0)) <
      ::pass::machine_precision) {
    // Uses quaternions
    const arma::vec& quaternion = {
        std::cos(roll_angle / 2.0) * std::cos(pitch_angle / 2.0) *
                std::cos(yaw_angle / 2.0) +
            std::sin(roll_angle / 2.0) * std::sin(pitch_angle / 2.0) *
                std::sin(yaw_angle / 2.0),
        std::sin(roll_angle / 2.0) * std::cos(pitch_angle / 2.0) *
                std::cos(yaw_angle / 2.0) -
            std::cos(roll_angle / 2.0) * std::sin(pitch_angle / 2.0) *
                std::sin(yaw_angle / 2.0),
        std::cos(roll_angle / 2.0) * std::sin(pitch_angle / 2.0) *
                std::cos(yaw_angle / 2.0) +
            std::sin(roll_angle / 2.0) * std::cos(pitch_angle / 2.0) *
                std::sin(yaw_angle / 2.0),
        std::cos(roll_angle / 2.0) * std::cos(pitch_angle / 2.0) *
                std::sin(yaw_angle / 2.0) -
            std::sin(roll_angle / 2.0) * std::sin(pitch_angle / 2.0) *
                std::cos(yaw_angle / 2.0)};
    return arma::mat(
        {{1.0 - 2.0 * (std::pow(quaternion(2), 2.0) +
                       std::pow(quaternion(3), 2.0)),
          2.0 * (quaternion(1) * quaternion(2) - quaternion(3) * quaternion(0)),
          2.0 *
              (quaternion(1) * quaternion(3) + quaternion(2) * quaternion(0))},
         {2.0 * (quaternion(1) * quaternion(2) + quaternion(3) * quaternion(0)),
          1.0 - 2.0 * (std::pow(quaternion(1), 2.0) +
                       std::pow(quaternion(3), 2.0)),
          2.0 *
              (quaternion(2) * quaternion(3) - quaternion(1) * quaternion(0))},
         {2.0 * (quaternion(1) * quaternion(3) - quaternion(2) * quaternion(0)),
          2.0 * (quaternion(2) * quaternion(3) + quaternion(1) * quaternion(0)),
          1.0 - 2.0 * (std::pow(quaternion(1), 2.0) +
                       std::pow(quaternion(2), 2.0))}});
  } else {
    // Uses Z,Y,X Tait-Bryan angles
    return arma::mat(
        {{std::cos(yaw_angle) * std::cos(pitch_angle),
          std::cos(yaw_angle) * std::sin(pitch_angle) * std::sin(roll_angle) -
              std::sin(yaw_angle) * std::cos(roll_angle),
          std::cos(yaw_angle) * std::sin(pitch_angle) * std::cos(roll_angle) +
              std::sin(yaw_angle) * std::sin(roll_angle)},
         {std::sin(yaw_angle) * std::cos(pitch_angle),
          std::sin(yaw_angle) * std::sin(pitch_angle) * std::sin(roll_angle) +
              std::cos(yaw_angle) * std::cos(roll_angle),
          std::sin(yaw_angle) * std::sin(pitch_angle) * std::cos(roll_angle) -
              std::cos(yaw_angle) * std::sin(roll_angle)},
         {-std::sin(pitch_angle), std::cos(pitch_angle) * std::sin(roll_angle),
          std::cos(pitch_angle) * std::cos(roll_angle)}});
  }
}

std::vector<arma::vec::fixed<2>> circle_circle_intersections(
    const arma::vec::fixed<2>& first_centre, const double first_radius,
    const arma::vec::fixed<2>& second_centre, const double second_radius) {
  assert(first_centre.is_finite() &&
         "circle_circle_intersections: The first centre must be finite.");
  assert(first_radius >= 0 &&
         "circle_circle_intersections: The first radius must be positive "
         "(including 0).");
  assert(second_centre.is_finite() &&
         "circle_circle_intersections: The second centre must be finite.");
  assert(second_radius >= 0 &&
         "circle_circle_intersections: The second radius must be positive "
         "(including 0).");

  /* The circle circle intersection points are calculated as follows:
   *  0. We assume that both centers are on the x-axis and `first_centre` is at
   *     (0, 0).
   *     **Note:** This assumptions are lifted later on.
   *  1. Calculate the distance between both centers.
   *  2. Calculate the x-part of the intersection point.
   *
   *   y-axis
   *      ^
   *      |
   *      |
   *      |----- distance ----|
   *      |                   |
   *      +-----------------------------> x-axis
   *      ^         ^         ^
   *      |         |         |
   * `first_center` |   `second_centrer`
   *                |
   *        intersection point (only the x-part)
   *
   * After the first to steps we should get a picture like above. At this point,
   * we either got 2 intersection points, when the distance between the x-part
   * of the intersection point and `first_center` is less than `first_radius`, 1
   * intersection point if both are equal and no intersection if its greater.
   *
   * A: We got 1 intersection
   *   3. Generate a unit vector pointing from `first_centre` to
   *      `second_centre`. 4. Set the y-part to 0. 5. Scale, rotate and
   *      translate the (x, y)-coordinate to be within the actual coordinate
   *      system (removing the assumption from 1.)
   *
   * B: We got 2 intersections
   *   3. Generate a unit vector pointing from `first_centre` to
   *      `second_centre`. 4. Calculate the y-part. 5. Scale, rotate and
   *      translate the (x, y)-coordinate to be within the actual coordinate
   *      system (removing the assumption from 1.)
   */

  const double distance = arma::norm(second_centre - first_centre);
  if (distance < ::pass::machine_precision &&
      std::abs(first_radius - second_radius) < ::pass::machine_precision) {
    // Both circles are identical ...
    if (first_radius > ::pass::machine_precision) {
      // ... but have a non-zero radius.
      throw std::range_error(
          "circle_circle_intersections: Both centers and radii (> 0) are "
          "identical, resulting in infinite intersections.");
    }
    // ... and dots.
    return {first_centre};
  }

  if (distance > first_radius + second_radius ||
      distance < std::abs(first_radius - second_radius)) {
    // Both circles are either too far away or too close.
    return {};
  }

  const double x = (std::pow(first_radius, 2.0) - std::pow(second_radius, 2.0) +
                    std::pow(distance, 2.0)) /
                   (2.0 * distance);
  const arma::vec::fixed<2>& unit_vector =
      (second_centre - first_centre) / distance;

  if (std::abs(first_radius - std::abs(x)) < ::pass::machine_precision) {
    // One intersection
    return std::vector<arma::vec::fixed<2>>({{first_centre + unit_vector * x}});
  } else {
    // Two intersections
    const double y = std::sqrt(std::pow(first_radius, 2.0) - std::pow(x, 2.0));
    return std::vector<arma::vec::fixed<2>>(
        {{first_centre(0) + unit_vector(0) * x - unit_vector(1) * y,
          first_centre(1) + unit_vector(1) * x + unit_vector(0) * y},
         {first_centre(0) + unit_vector(0) * x + unit_vector(1) * y,
          first_centre(1) + unit_vector(1) * x - unit_vector(0) * y}});
  }
}

std::vector<arma::vec::fixed<3>> circle_sphere_intersections(
    const arma::vec::fixed<3>& circle_centre, const double circle_radius,
    const arma::vec::fixed<3>& circle_normal,
    const arma::vec::fixed<3>& sphere_centre, const double sphere_radius) {
  assert(circle_centre.is_finite() &&
         "circle_sphere_intersections: The circle centre must be finite.");
  assert(circle_radius >= 0 &&
         "circle_sphere_intersections: The circle radius must be positive "
         "(including 0).");
  assert(circle_normal.is_finite() &&
         "circle_sphere_intersections: The circle normal must be finite.");
  assert(sphere_centre.is_finite() &&
         "circle_sphere_intersections: The sphere centre must be finite.");
  assert(sphere_radius >= 0 &&
         "circle_sphere_intersections: The sphere radius must be positive "
         "(including 0).");

  /* The circle sphere intersection points are calculated as follows:
   * 0. We assume that both centers are on the x-axis, ´circle_normal` is
   * perpendicular to the x- and y-axis and `circle_centre` is at (0, 0, 0).
   *    **Note:** This assumptions are lifted later on.
   * 1. Calculate the shortest distance between the sphere's centre and circle's
   * plane. 2. Calculate the centre of the sphere's inner circle segment, placed
   * an the same plane as the given circle. 3. Calculate the radius of the
   * sphere's inner circle segment. Given the provided circle and the sphere's
   * inner circle, the problem is now reduced to a circle circle
   * intersection. 4. Calculate the x-part of the intersection point. At this
   * point, we either got 2 intersection points, when the distance between
   * the x-part of the intersection point and `circle_centre` is less than
   * `circle_radius`, 1 intersection point if both are equal and no intersection
   * if its greater.
   *
   * A: We got 1 intersection
   *   5. Generate a unit vector pointing from `circle_centre` to
   * `inner_centre`. 6. Set the y-part to 0. 7. Scale and translate the unit
   * vector to be within the actual coordinate system (removing the assumption
   * from 1.).
   *
   * B: We got 2 intersections
   *   5. Generate a unit vector pointing from `circle_centre` to
   * `inner_centre`. 6. Calculate the y-part. 7. Generate a second unit vector,
   * perpendicular to `x_unit_vector` and `circle_normal`. 8. Scale and
   * translate both unit vectors to be within the actual coordinate system
   * (removing the assumption from 1.).
   */

  const double inner_distance =
      arma::dot(circle_normal, circle_centre - sphere_centre);
  if (std::abs(inner_distance) - sphere_radius > ::pass::machine_precision) {
    // The circle's plane does not intersect with the sphere.
    return {};
  };

  const arma::vec::fixed<3>& inner_centre =
      sphere_centre + inner_distance * circle_normal;
  // Due to rounding errors, the inner result might be negative instead of being
  // 0.
  const double inner_radius = std::sqrt(std::max(
      0.0, std::pow(sphere_radius, 2.0) - std::pow(inner_distance, 2.0)));

  const double distance = arma::norm(inner_centre - circle_centre);
  if (distance < ::pass::machine_precision &&
      std::abs(inner_radius - circle_radius) < ::pass::machine_precision) {
    // Both circles are identical ...
    if (circle_radius > ::pass::machine_precision) {
      // ... but have a non-zero radius.
      throw std::range_error(
          "circle_sphere_intersections: Both centers and radii (> 0) are "
          "identical, resulting in infinite intersections.");
    }

    // ... and dots.
    return {circle_centre};
  }

  if (distance - circle_radius - inner_radius > ::pass::machine_precision ||
      std::abs(circle_radius - inner_radius) - distance >=
          ::pass::machine_precision) {
    // Both circles are either to far away or to close.
    return {};
  };

  const double x = (std::pow(circle_radius, 2.0) - std::pow(inner_radius, 2.0) +
                    std::pow(distance, 2.0)) /
                   (2.0 * distance);
  const arma::vec::fixed<3>& x_unit_vector =
      (inner_centre - circle_centre) / distance;

  if (std::abs(circle_radius - std::abs(x)) < ::pass::machine_precision) {
    // One intersection
    return std::vector<arma::vec::fixed<3>>(
        {circle_centre + x * x_unit_vector});
  } else {
    // Two intersections
    const double y = std::sqrt(std::pow(circle_radius, 2.0) - std::pow(x, 2.0));
    const arma::vec::fixed<3>& y_unit_vector =
        arma::normalise(arma::cross(x_unit_vector, circle_normal));
    return std::vector<arma::vec::fixed<3>>(
        {circle_centre + x * x_unit_vector + y * y_unit_vector,
         circle_centre + x * x_unit_vector - y * y_unit_vector});
  }
}
}  // namespace pass
