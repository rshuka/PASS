// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's            //
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //

#include "pass_bits/helper/astro_problems/pow_swing_by_inv.hpp"
#include <cmath>

void pass::pow_swing_by_inv(const double Vin, const double Vout, const double alpha,
                            double &DV, double &rp)
{
  const int maxiter = 30;
  int i = 0;
  double err = 1.0;
  const double tolerance = 1e-8;

  double aIN = 1.0 / pow(Vin, 2);   // semimajor axis of the incoming hyperbola
  double aOUT = 1.0 / pow(Vout, 2); // semimajor axis of the incoming hyperbola

  rp = 1.0;
  while ((err > tolerance) && (i < maxiter))
  {
    i++;
    double f = asin(aIN / (aIN + rp)) + asin(aOUT / (aOUT + rp)) - alpha;
    double df = -aIN / sqrt((rp + 2 * aIN) * rp) / (aIN + rp) -
                aOUT / sqrt((rp + 2 * aOUT) * rp) / (aOUT + rp);
    double rp_new = rp - f / df;
    if (rp_new > 0)
    {
      err = fabs(rp_new - rp);
      rp = rp_new;
    }
    else
      rp /= 2.0;
  }

  // Evaluation of the DV
  DV = fabs(sqrt(Vout * Vout + (2.0 / rp)) - sqrt(Vin * Vin + (2.0 / rp)));
}
