#include "pass_bits/helper/astro_problems/constants.hpp"
#include "pass_bits/helper/astro_problems/astro_helpers.hpp"
#include <cmath>
#include <iterator>
#include <numeric>

namespace pass
{
const celestial_body celestial_body::SUN(1.32712428e11, 0, 0);
const celestial_body celestial_body::MERCURY(22321, 0, 0);
const celestial_body celestial_body::VENUS(324860, 6351.8, 0.01);
const celestial_body celestial_body::EARTH(398601.19, 6778.1, 0.01);
const celestial_body celestial_body::MARS(42828.3, 6000, 0.01);
const celestial_body celestial_body::JUPITER(126.7e6, 600000, 0.001);
const celestial_body celestial_body::SATURN(37.9e6, 70000, 0.01);
const celestial_body celestial_body::URANUS(5.78e6, 0, 0);
const celestial_body celestial_body::NEPTUNE(6.8e6, 0, 0);

celestial_body::celestial_body(double mu, double penalty,
                               double penalty_coefficient) noexcept
    : mu(mu), penalty(penalty), penalty_coefficient(penalty_coefficient) {}

std::pair<std::array<double, 3>, std::array<double, 3>>
celestial_body::ephemeris(const double mjd2000) const
{
  const double pi = acos(-1.0);
  const double RAD = pi / 180.0;
  const double AU = 149597870.66; // Astronomical Unit
  const double KM = AU;
  std::array<double, 6> Kepl_Par;
  double XM;

  double T = (mjd2000 + 36525.00) / 36525.00;

  if (this == &celestial_body::MERCURY)
  {
    Kepl_Par[0] = (0.38709860);
    Kepl_Par[1] = (0.205614210 + 0.000020460 * T - 0.000000030 * T * T);
    Kepl_Par[2] = (7.002880555555555560 + 1.86083333333333333e-3 * T -
                   1.83333333333333333e-5 * T * T);
    Kepl_Par[3] = (4.71459444444444444e+1 + 1.185208333333333330 * T +
                   1.73888888888888889e-4 * T * T);
    Kepl_Par[4] = (2.87537527777777778e+1 + 3.70280555555555556e-1 * T +
                   1.20833333333333333e-4 * T * T);
    XM = 1.49472515288888889e+5 + 6.38888888888888889e-6 * T;
    Kepl_Par[5] = (1.02279380555555556e2 + XM * T);
  }
  else if (this == &celestial_body::VENUS)
  {
    Kepl_Par[0] = (0.72333160);
    Kepl_Par[1] = (0.006820690 - 0.000047740 * T + 0.0000000910 * T * T);
    Kepl_Par[2] = (3.393630555555555560 + 1.00583333333333333e-3 * T -
                   9.72222222222222222e-7 * T * T);
    Kepl_Par[3] = (7.57796472222222222e+1 + 8.9985e-1 * T + 4.1e-4 * T * T);
    Kepl_Par[4] = (5.43841861111111111e+1 + 5.08186111111111111e-1 * T -
                   1.38638888888888889e-3 * T * T);
    XM = 5.8517803875e+4 + 1.28605555555555556e-3 * T;
    Kepl_Par[5] = (2.12603219444444444e2 + XM * T);
  }
  else if (this == &celestial_body::EARTH)
  {
    Kepl_Par[0] = (1.000000230);
    Kepl_Par[1] = (0.016751040 - 0.000041800 * T - 0.0000001260 * T * T);
    Kepl_Par[2] = (0.00);
    Kepl_Par[3] = (0.00);
    Kepl_Par[4] =
        (1.01220833333333333e+2 + 1.7191750 * T +
         4.52777777777777778e-4 * T * T + 3.33333333333333333e-6 * T * T * T);
    XM = 3.599904975e+4 - 1.50277777777777778e-4 * T -
         3.33333333333333333e-6 * T * T;
    Kepl_Par[5] = (3.58475844444444444e2 + XM * T);
  }
  else if (this == &celestial_body::MARS)
  {
    Kepl_Par[0] = (1.5236883990);
    Kepl_Par[1] = (0.093312900 + 0.0000920640 * T - 0.0000000770 * T * T);
    Kepl_Par[2] =
        (1.850333333333333330 - 6.75e-4 * T + 1.26111111111111111e-5 * T * T);
    Kepl_Par[3] =
        (4.87864416666666667e+1 + 7.70991666666666667e-1 * T -
         1.38888888888888889e-6 * T * T - 5.33333333333333333e-6 * T * T * T);
    Kepl_Par[4] = (2.85431761111111111e+2 + 1.069766666666666670 * T +
                   1.3125e-4 * T * T + 4.13888888888888889e-6 * T * T * T);
    XM = 1.91398585e+4 + 1.80805555555555556e-4 * T +
         1.19444444444444444e-6 * T * T;
    Kepl_Par[5] = (3.19529425e2 + XM * T);
  }
  else if (this == &celestial_body::JUPITER)
  {
    Kepl_Par[0] = (5.2025610);
    Kepl_Par[1] = (0.048334750 + 0.000164180 * T - 0.00000046760 * T * T -
                   0.00000000170 * T * T * T);
    Kepl_Par[2] = (1.308736111111111110 - 5.69611111111111111e-3 * T +
                   3.88888888888888889e-6 * T * T);
    Kepl_Par[3] =
        (9.94433861111111111e+1 + 1.010530 * T +
         3.52222222222222222e-4 * T * T - 8.51111111111111111e-6 * T * T * T);
    Kepl_Par[4] = (2.73277541666666667e+2 + 5.99431666666666667e-1 * T +
                   7.0405e-4 * T * T + 5.07777777777777778e-6 * T * T * T);
    XM = 3.03469202388888889e+3 - 7.21588888888888889e-4 * T +
         1.78444444444444444e-6 * T * T;
    Kepl_Par[5] = (2.25328327777777778e2 + XM * T);
  }
  else if (this == &celestial_body::SATURN)
  {
    Kepl_Par[0] = (9.5547470);
    Kepl_Par[1] = (0.055892320 - 0.00034550 * T - 0.0000007280 * T * T +
                   0.000000000740 * T * T * T);
    Kepl_Par[2] =
        (2.492519444444444440 - 3.91888888888888889e-3 * T -
         1.54888888888888889e-5 * T * T + 4.44444444444444444e-8 * T * T * T);
    Kepl_Par[3] =
        (1.12790388888888889e+2 + 8.73195138888888889e-1 * T -
         1.52180555555555556e-4 * T * T - 5.30555555555555556e-6 * T * T * T);
    Kepl_Par[4] =
        (3.38307772222222222e+2 + 1.085220694444444440 * T +
         9.78541666666666667e-4 * T * T + 9.91666666666666667e-6 * T * T * T);
    XM = 1.22155146777777778e+3 - 5.01819444444444444e-4 * T -
         5.19444444444444444e-6 * T * T;
    Kepl_Par[5] = (1.75466216666666667e2 + XM * T);
  }
  else if (this == &celestial_body::URANUS)
  {
    Kepl_Par[0] = (19.218140);
    Kepl_Par[1] = (0.04634440 - 0.000026580 * T + 0.0000000770 * T * T);
    Kepl_Par[2] =
        (7.72463888888888889e-1 + 6.25277777777777778e-4 * T + 3.95e-5 * T * T);
    Kepl_Par[3] = (7.34770972222222222e+1 + 4.98667777777777778e-1 * T +
                   1.31166666666666667e-3 * T * T);
    Kepl_Par[4] =
        (9.80715527777777778e+1 + 9.85765e-1 * T -
         1.07447222222222222e-3 * T * T - 6.05555555555555556e-7 * T * T * T);
    XM = 4.28379113055555556e+2 + 7.88444444444444444e-5 * T +
         1.11111111111111111e-9 * T * T;
    Kepl_Par[5] = (7.26488194444444444e1 + XM * T);
  }
  else if (this == &celestial_body::NEPTUNE)
  {
    Kepl_Par[0] = (30.109570);
    Kepl_Par[1] = (0.008997040 + 0.0000063300 * T - 0.0000000020 * T * T);
    Kepl_Par[2] = (1.779241666666666670 - 9.54361111111111111e-3 * T -
                   9.11111111111111111e-6 * T * T);
    Kepl_Par[3] =
        (1.30681358333333333e+2 + 1.0989350 * T +
         2.49866666666666667e-4 * T * T - 4.71777777777777778e-6 * T * T * T);
    Kepl_Par[4] = (2.76045966666666667e+2 + 3.25639444444444444e-1 * T +
                   1.4095e-4 * T * T + 4.11333333333333333e-6 * T * T * T);
    XM = 2.18461339722222222e+2 - 7.03333333333333333e-5 * T;
    Kepl_Par[5] = (3.77306694444444444e1 + XM * T);
  }
  else
  {
    throw std::domain_error(
        "`celestial_body::conversion` not implemented for SUN");
  }

  // conversion of AU into KM
  Kepl_Par[0] *= KM;

  // conversion of DEG into RAD
  Kepl_Par[2] *= RAD;
  Kepl_Par[3] *= RAD;
  Kepl_Par[4] *= RAD;
  Kepl_Par[5] *= RAD;
  Kepl_Par[5] = fmod(Kepl_Par[5], 2.0 * pi);

  // Conversion from Mean Anomaly to Eccentric Anomaly via Kepler's equation
  Kepl_Par[5] = mean_to_eccentric(Kepl_Par[5], Kepl_Par[1]);

  // Position and Velocity evaluation according to j2000 system
  return conversion(Kepl_Par, celestial_body::SUN.mu);
}

asteroid::asteroid(std::array<double, 6> keplerian, double epoch,
                   double mu) noexcept
    : keplerian(keplerian), epoch(epoch), mu(mu) {}

std::pair<std::array<double, 3>, std::array<double, 3>> asteroid::ephemeris(
    const double jd) const
{
  const double pi = acos(-1.0);
  const double RAD = pi / 180.0;
  const double AU = 149597870.66; // Astronomical Unit
  double a, e, i, W, w, M, jdepoch, DT, n, E;
  std::array<double, 6> V;

  a = keplerian[0] * AU; // in km
  e = keplerian[1];
  i = keplerian[2];
  W = keplerian[3];
  w = keplerian[4];
  M = keplerian[5];
  jdepoch = epoch + 2400000.5;
  DT = (jd - jdepoch) * 86400;
  n = sqrt(celestial_body::SUN.mu / pow(a, 3));

  M = M / 180.0 * pi;
  M += n * DT;
  M = fmod(M, 2 * pi);
  E = mean_to_eccentric(M, e);
  V[0] = a;
  V[1] = e;
  V[2] = i * RAD;
  V[3] = W * RAD;
  V[4] = w * RAD;
  V[5] = E;

  return conversion(V, celestial_body::SUN.mu);
}
} // namespace pass
