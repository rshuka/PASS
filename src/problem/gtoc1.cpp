#include "pass_bits/problem/gtoc1.hpp"

// assert
#include <cassert>

#include "pass_bits/helper/gtoc1/astro_helpers.hpp"
#include "pass_bits/helper/gtoc1/vector3d_helpers.hpp"

pass::gtoc1::gtoc1()
    : sequence({&celestial_body::EARTH, &celestial_body::VENUS,
                &celestial_body::EARTH, &celestial_body::VENUS,
                &celestial_body::EARTH, &celestial_body::JUPITER,
                &celestial_body::SATURN}),
      rev_flag({0, 0, 0, 0, 0, 0, 1, 0}),
      destination(
          {{2.5897261, 0.2734625, 6.40734, 128.34711, 264.78691, 320.479555},
           53600.0,
           0.0}),
      Isp(2500.0),
      mass(1500.0),
      DVlaunch(2.5),
      problem({3000, 14, 14, 14, 14, 100, 366, 300},
              {10000, 2000, 2000, 2000, 2000, 9000, 9000, 9000}) {}

double pass::gtoc1::evaluate(const arma::vec& parameter) const {
  assert(parameter.n_elem == dimension() &&
         "`parameter` has incompatible dimension");

  std::array<double, 6> rp{};
  std::array<double, 8> DV{};
  const int n = 8;

  const double g = 9.80665 / 1000.0;  // Gravity

  // r and  v in heliocentric coordinate system
  // position
  std::array<std::array<double, 3>, 8> r;
  // velocity
  std::array<std::array<double, 3>, 8> v;

  {
    double totalTime = 0;
    for (size_t i = 0; i < 7; i++) {
      totalTime += parameter[i];
      auto result = sequence[i]->ephemeris(totalTime);
      r[i] = result.first;
      v[i] = result.second;
    }
    totalTime += parameter[7];
    auto result = destination.ephemeris(totalTime + 2451544.5);
    r[7] = result.first;
    v[7] = result.second;
  }

  std::array<double, 3> current_section_departure_velocity;
  std::array<double, 3> current_section_arrival_velocity;
  for (size_t i = 0; i <= n - 2; i++) {
    std::array<double, 3> last_section_departure_velocity =
        current_section_departure_velocity;
    std::array<double, 3> last_section_arrival_velocity =
        current_section_arrival_velocity;

    bool longWay = pass::gtoc::cross_product(r[i], r[i + 1])[2] > 0
                       ? rev_flag[i]
                       : !rev_flag[i];

    const auto lambert_solution =
        lambert(r[i], r[i + 1], parameter[i + 1] * 24 * 60 * 60,
                celestial_body::SUN.mu, longWay);
    current_section_departure_velocity = lambert_solution.departure_velocity;
    current_section_arrival_velocity = lambert_solution.arrival_velocity;

    if (i == 0) {
      // Earth launch
      DV[0] = pass::gtoc::norm(
          pass::gtoc::sub(current_section_departure_velocity, v[0]));
    } else {
      double Vin = pass::gtoc::norm(
          pass::gtoc::sub(last_section_arrival_velocity, v[i]));
      double Vout = pass::gtoc::norm(
          pass::gtoc::sub(current_section_departure_velocity, v[i]));

      // calculation of delta V at pericenter
      const auto swing_by_solution = PowSwingByInv(
          Vin, Vout,
          acos(pass::gtoc::dot_product(
                   pass::gtoc::sub(last_section_arrival_velocity, v[i]),
                   pass::gtoc::sub(current_section_departure_velocity, v[i])) /
               (Vin * Vout)));
      DV[i] = swing_by_solution.first;
      rp[i - 1] = swing_by_solution.second;
      rp[i - 1] *= sequence[i]->mu;
    }
  }

  double DVtot = std::accumulate<double*, double>(std::next(DV.begin()),
                                                  std::prev(DV.end()), 0);

  // Build Penalty
  for (size_t i = 0; i < n - 2; i++) {
    if (rp[i] < sequence[i + 1]->penalty)
      DVtot += sequence[i + 1]->penalty_coefficient *
               fabs(rp[i] - sequence[i + 1]->penalty);
  }

  // Launcher Constraint
  if (DV[0] > DVlaunch) DVtot += (DV[0] - DVlaunch);

  // Evaluation of satellite final mass
  double final_mass = mass * exp(-DVtot / (Isp * g));

  // arrival relative velocity at the asteroid;
  auto arrival_velocity =
      pass::gtoc::sub(v[n - 1], current_section_arrival_velocity);

  return -final_mass *
         fabs(pass::gtoc::dot_product(arrival_velocity, v[n - 1]));
}
