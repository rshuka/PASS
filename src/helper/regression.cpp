#include "pass_bits/helper/regression.hpp"
#include <math.h> // sqrt, pow

arma::colvec pass::regression::linear_model(const int elements, const arma::colvec x_values, const arma::colvec y_values)
{

  arma::colvec model(3);

  double sum_x = 0.0;  /* sum of x     */
  double sum_x2 = 0.0; /* sum of x^2  */
  double sum_xy = 0.0; /* sum of x * y */
  double sum_y = 0.0;  /* sum of y     */
  double sum_y2 = 0.0; /* sum of y^2  */

  for (int i = 0; i < elements; i++)
  {
    sum_x += x_values(i);
    sum_x2 += std::pow(x_values(i), 2);
    sum_xy += x_values(i) * y_values(i);
    sum_y += y_values(i);
    sum_y2 += std::pow(y_values(i), 2);
  }

  double denom = elements * sum_x2 - std::pow(sum_x, 2);

  if (denom == 0)
  {
    model(0) = 0;
    model(1) = 0;
    model(2) = 0;
    throw std::runtime_error(
        "You have a singular matrix and the problem can't be solved");
  }

  model(0) = (elements * sum_xy - sum_x * sum_y) / denom;
  model(1) = (sum_y * sum_x2 - sum_x * sum_xy) / denom;

  model(2) = (sum_xy - sum_x * sum_y / elements) / /* compute correlation coeff */
             std::sqrt((sum_x2 - std::pow(sum_x, 2) / elements) *
                       (sum_y2 - std::pow(sum_y, 2) / elements));

  return model;
}

double pass::regression::predict_linear(const double x, const arma::colvec model)
{
  return model(0) * x + model(1);
}
