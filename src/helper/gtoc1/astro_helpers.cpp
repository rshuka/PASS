#include "pass_bits/helper/gtoc1/astro_helpers.hpp"

#include <cmath>
#include <stdexcept>

#include "pass_bits/helper/gtoc1/Astro_Functions.hpp"
#include "pass_bits/helper/gtoc1/vector3d_helpers.hpp"

pass::lambert_solution pass::lambert(std::array<double, 3> r1,
                                     std::array<double, 3> r2, double t,
                                     const double mu, const int lw) {
  if (t <= 0) {
    throw std::invalid_argument(
        "ERROR in Lambert Solver: Negative Time in input.");
  }

  double R = sqrt(pass::gtoc::dot_product(r1, r1));
  double V = sqrt(mu / R);
  double T = R / V;

  // working with non-dimensional radii and time-of-flight
  t /= T;
  r1 = pass::gtoc::div(r1, R);
  r2 = pass::gtoc::div(r2, R);

  // Evaluation of the relevant geometry parameters in non dimensional units
  // R2 module
  double r2_mod = sqrt(pass::gtoc::dot_product(r2, r2));

  double theta = acos(pass::gtoc::dot_product(r1, r2) / r2_mod);

  if (lw) theta = 2 * acos(-1.0) - theta;

  // non-dimensional chord
  double c = sqrt(1 + r2_mod * (r2_mod - 2.0 * cos(theta)));

  // non dimesnional semi-perimeter
  double s = (1 + r2_mod + c) / 2.0;

  // minimum energy ellipse semi major axis
  double am = s / 2.0;

  // lambda parameter defined in Battin's Book
  double lambda = sqrt(r2_mod) * cos(theta / 2.0) / s;

  double x;
  {
    // We start finding the log(x+1) value of the solution conic:
    // NO MULTI REV --> (1 SOL)
    //	inn1=-.5233;    //first guess point
    //  inn2=.5233;     //second guess point
    double x1 = log(0.4767);
    double x2 = log(1.5233);
    double y1 = log(x2tof(-.5233, s, c, lw)) - log(t);
    double y2 = log(x2tof(.5233, s, c, lw)) - log(t);

    // Increasing the tolerance does not bring any advantage as the
    // precision is usually greater anyway (due to the rectification of the tof
    // graph) except near particular cases such as parabolas in which cases a
    // lower precision allow for usual convergence.
    const double tolerance = 1e-11;
    double error = 1;

    // Newton iterations
    double x_new = 0.0, y_new;
    while ((error > tolerance) && (y1 != y2)) {
      x_new = (x1 * y2 - y1 * x2) / (y2 - y1);
      y_new = logf(x2tof(expf(x_new) - 1, s, c, lw)) - logf(t);
      x1 = x2;
      y1 = y2;
      x2 = x_new;
      y2 = y_new;
      error = fabs(x1 - x_new);
    }
    x = expf(x_new) - 1;
  }

  // The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
  // now need the conic. As for transfer angles near to pi the lagrange
  // coefficient technique goes singular (dg approaches a zero/zero that is
  // numerically bad) we here use a different technique for those cases. When
  // the transfer angle is exactly equal to pi, then the ih unit vector is not
  // determined. The remaining equations, though, are still valid.

  double a = am / (1 - x * x);

  double eta, eta_squared;
  // psi evaluation
  if (x < 1) {
    // ellipse
    double beta = 2 * asin(sqrt((s - c) / (2 * a)));
    if (lw) beta = -beta;
    double alfa = 2 * acos(x);
    double psi = (alfa - beta) / 2;
    eta_squared = 2 * a * pow(sin(psi), 2) / s;
    eta = sqrt(eta_squared);
  } else {
    // hyperbola
    double beta = 2 * asinh(sqrt((c - s) / (2 * a)));
    if (lw) beta = -beta;
    double alfa = 2 * acosh(x);
    double psi = (alfa - beta) / 2;
    eta_squared = -2 * a * pow(sinh(psi), 2) / s;
    eta = sqrt(eta_squared);
  }

  // parameter of the solution
  double p = (r2_mod / (am * eta_squared)) * pow(sin(theta / 2), 2);

  std::array<double, 3> ih =
      pass::gtoc::unit_vector(pass::gtoc::cross_product(r1, r2));
  if (lw) ih = pass::gtoc::mul(ih, -1);

  pass::lambert_solution solution;

  double vr1 = (1 / (eta * sqrt(am))) * (2 * lambda * am - (lambda + x * eta));
  double vt1 = sqrt(p);
  solution.departure_velocity = pass::gtoc::mul(
      pass::gtoc::add(pass::gtoc::mul(r1, vr1),
                      pass::gtoc::mul(pass::gtoc::cross_product(ih, r1), vt1)),
      V);

  double vt2 = vt1 / r2_mod;
  double vr2 = -vr1 + (vt1 - vt2) / tan(theta / 2);
  solution.arrival_velocity = pass::gtoc::mul(
      pass::gtoc::add(
          pass::gtoc::div(pass::gtoc::mul(r2, vr2), r2_mod),
          pass::gtoc::mul(
              pass::gtoc::cross_product(ih, pass::gtoc::unit_vector(r2)), vt2)),
      V);

  return solution;
}

std::pair<double, double> pass::PowSwingByInv(const double Vin,
                                              const double Vout,
                                              const double alpha) {
  double DV, rp;
  const int maxiter = 30;
  int i = 0;
  double err = 1.0;
  double f, df;  // function and its derivative
  double rp_new;
  const double tolerance = 1e-8;

  double aIN = 1.0 / pow(Vin, 2);    // semimajor axis of the incoming hyperbola
  double aOUT = 1.0 / pow(Vout, 2);  // semimajor axis of the incoming hyperbola

  rp = 1.0;
  while ((err > tolerance) && (i < maxiter)) {
    i++;
    f = asin(aIN / (aIN + rp)) + asin(aOUT / (aOUT + rp)) - alpha;
    df = -aIN / sqrt((rp + 2 * aIN) * rp) / (aIN + rp) -
         aOUT / sqrt((rp + 2 * aOUT) * rp) / (aOUT + rp);
    rp_new = rp - f / df;
    if (rp_new > 0) {
      err = fabs(rp_new - rp);
      rp = rp_new;
    } else
      rp /= 2.0;
  }

  // Evaluation of the DV
  DV = fabs(sqrt(Vout * Vout + (2.0 / rp)) - sqrt(Vin * Vin + (2.0 / rp)));

  return {DV, rp};
}

std::pair<std::array<double, 3>, std::array<double, 3>> pass::conversion(
    const std::array<double, 6>& E, const double mu) {
  double a, e, i, omg, omp, theta;
  double b, n;
  double X_per[3], X_dotper[3];
  double R[3][3];

  a = E[0];
  e = E[1];
  i = E[2];
  omg = E[3];
  omp = E[4];
  theta = E[5];

  b = a * sqrt(1 - e * e);
  n = sqrt(mu / pow(a, 3));

  const double sin_theta = sin(theta);
  const double cos_theta = cos(theta);

  X_per[0] = a * (cos_theta - e);
  X_per[1] = b * sin_theta;

  X_dotper[0] = -(a * n * sin_theta) / (1 - e * cos_theta);
  X_dotper[1] = (b * n * cos_theta) / (1 - e * cos_theta);

  const double cosomg = cos(omg);
  const double cosomp = cos(omp);
  const double sinomg = sin(omg);
  const double sinomp = sin(omp);
  const double cosi = cos(i);
  const double sini = sin(i);

  R[0][0] = cosomg * cosomp - sinomg * sinomp * cosi;
  R[0][1] = -cosomg * sinomp - sinomg * cosomp * cosi;

  R[1][0] = sinomg * cosomp + cosomg * sinomp * cosi;
  R[1][1] = -sinomg * sinomp + cosomg * cosomp * cosi;

  R[2][0] = sinomp * sini;
  R[2][1] = cosomp * sini;

  std::pair<std::array<double, 3>, std::array<double, 3>> result;
  // evaluate position and velocity
  for (int i = 0; i < 3; i++) {
    result.first[i] = 0;
    result.second[i] = 0;
    for (int j = 0; j < 2; j++) {
      result.first[i] += R[i][j] * X_per[j];
      result.second[i] += R[i][j] * X_dotper[j];
    }
  }
  return result;
}

double pass::mean_to_eccentric(const double m, const double e) {
  return Mean2Eccentric(m, e);
}
