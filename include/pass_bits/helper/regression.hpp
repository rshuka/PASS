#pragma once
#include <array>
#include <armadillo>
#include <cassert>
namespace pass
{
/**
   * Regression Implementation
   *
   * Linear and polynomial regression
   */
class regression
{
public:
  /**
   *  Creates the linear model for the given data set
   *  using the linear least squares method
   *  https://en.wikipedia.org/wiki/Linear_least_squares
   *
   * @return a rowvec like y = mx + b and the correlation coeff r^2
   * first element -> m
   * second element -> b
   * thirt element -> r^2
   *
   * r^2 in range [0, 1]
   * -> 0 worst fit
   * -> 1 best fit
   */
  arma::rowvec linear_model(const arma::rowvec &x_values, const arma::rowvec &y_values);

  /**
   *  Creates the polynomial model (3 Grade!!) for the given data set
   *  using the least squares method
   *  https://en.wikipedia.org/wiki/Linear_least_squares
   *
   * Code based on Manas Sharma example
   * @see https://www.bragitoff.com/2015/09/c-program-for-polynomial-fit-least-squares/
   *
   * @return a rowvec like f(x) = p1*x^3 + p2*x^2 + p3*x + p4
   * first element -> p4
   * second element -> p4
   * thirt element -> p2
   * fourth element -> p1
   */

  arma::rowvec poly_model(const arma::rowvec &x_values, const arma::rowvec &y_values);

  /**
   * Predict a value x using the linear model
   */
  double predict_linear(const double &x, const arma::rowvec &model);

  /**
   * Predict a value x using the polynomial model
   */
  double predict_poly(const double &x, const arma::rowvec &model);
};
} // namespace pass
