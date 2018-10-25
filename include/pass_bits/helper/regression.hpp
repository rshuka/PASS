#pragma once
#include <armadillo>
namespace pass
{
/**
   * Regression Implementation
   *
   * Linear and polynomial regression
   *
   */
class regression
{
public:
  /**
   *  Creates the linear model for the given data set
   *  using the linear least squares method
   *  https://en.wikipedia.org/wiki/Linear_least_squares
   *
   * @return a colvec like y = mx + b and the correlation coeff r^2
   * first element -> m
   * second element -> b
   * thirt element -> r^2
   *
   * r^2 in range [0, 1]
   * -> 0 worst fit
   * -> 1 best fit
   */
  arma::colvec linear_model(const int elements, const arma::colvec x_values, const arma::colvec y_values);

  /**
   * Predict a value x using a linear model
   */

  double predict_linear(const double x, const arma::colvec model);
};
} // namespace pass
