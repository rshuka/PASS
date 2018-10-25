#include "pass_bits/helper/regression.hpp"
#include <math.h> // sqrt, pow

arma::rowvec pass::regression::linear_model(const arma::rowvec &x_values, const arma::rowvec &y_values)
{
  assert(x_values.n_elem > 1 && y_values.n_elem > 1 && "The training data should have more than one element");

  const int elements = x_values.n_elem;

  arma::rowvec model(3);

  double sum_x = 0.0;  // sum of x
  double sum_x2 = 0.0; // sum of x^2
  double sum_xy = 0.0; // sum of x * y
  double sum_y = 0.0;  // sum of y
  double sum_y2 = 0.0; // sum of y^2

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

  model(2) = (sum_xy - sum_x * sum_y / elements) /
             std::sqrt((sum_x2 - std::pow(sum_x, 2) / elements) *
                       (sum_y2 - std::pow(sum_y, 2) / elements));

  model(2) = std::pow(model(2), 2);

  return model;
}

arma::rowvec pass::regression::poly_model(const arma::rowvec &x_values, const arma::rowvec &y_values)
{
  assert(x_values.n_elem > 1 && y_values.n_elem > 1 && "The training data should have more than one element");

  int degree = 3; // degree of polynom

  arma::rowvec model(degree + 1);

  int i, j, k; // counters

  int elements = x_values.n_elem;

  //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
  arma::rowvec sigmaX(2 * degree + 1);

  //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
  arma::rowvec sigmaY(degree + 1);

  for (i = 0; i < 2 * degree + 1; i++)
  {
    sigmaX(i) = 0;

    for (j = 0; j < elements; j++)
    {
      //consecutive positions of the array will store Elements,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
      sigmaX(i) = sigmaX(i) + std::pow(x_values(j), i);
    }
  }

  //the Normal matrix(augmented) that will store the equations
  arma::mat matrix(degree + 1, degree + 2);

  //value of the final coefficients
  arma::rowvec coeff(degree + 1);

  for (i = 0; i <= degree; i++)
  {
    for (j = 0; j <= degree; j++)
    {
      //Build the Normal matrix by storing the corresponding coefficients at the right positions
      //except the last column of the matrix
      matrix(i, j) = sigmaX(i + j);
    }
  }

  for (i = 0; i < degree + 1; i++)
  {
    sigmaY(i) = 0;
    for (j = 0; j < elements; j++)
    {
      //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
      sigmaY(i) = sigmaY(i) + std::pow(x_values(j), i) * y_values(j);
    }
  }

  for (i = 0; i <= degree; i++)
  {
    //load the values of sigmaY as the last column of matrix(Normal Matrix but augmented)
    matrix(i, degree + 1) = sigmaY(i);
  }

  //'degree' is made 'degree+1' because the Gaussian Elimination part below was for 'degree' equations,
  //but here 'degree' is the degree of polynomial and for 'degree' degree we get 'degree'+1 equations
  degree = degree + 1;

  //From now Gaussian Elimination starts to solve the set of linear equations (Pivotisation)
  for (i = 0; i < degree; i++)
  {
    for (k = i + 1; k < degree; k++)
    {
      if (matrix(i, i) < matrix(k, i))
      {
        for (j = 0; j <= degree; j++)
        {
          double temp = matrix(i, j);
          matrix(i, j) = matrix(k, j);
          matrix(k, j) = temp;
        }
      }
    }
  }

  //loop to perform the gauss elimination
  for (i = 0; i < degree - 1; i++)
  {
    for (k = i + 1; k < degree; k++)
    {
      double t = matrix(k, i) / matrix(i, i);
      for (j = 0; j <= degree; j++)
      {
        //make the elements below the pivot elements equal to zero or elimnate the variables
        matrix(k, j) = matrix(k, j) - t * matrix(i, j);
      }
    }
  }

  for (i = degree - 1; i >= 0; i--) //back-substitution
  {
    //make the variable to be calculated equal to the rhs of the last equation
    coeff(i) = matrix(i, degree);

    for (j = 0; j < degree; j++)
    {
      //then subtract all the lhs values except the coefficient of the variable whose value is being calculated
      if (j != i)
      {
        coeff(i) = coeff(i) - matrix(i, j) * coeff(j);
      }
    }
    //now finally divide the rhs by the coefficient of the variable to be calculated
    coeff(i) = coeff(i) / matrix(i, i);
  }

  for (i = 0; i < degree; i++)
  {
    model(i) = coeff(i);
  }

  return model;
}

double pass::regression::predict_linear(const double &x, const arma::rowvec &model)
{
  assert(model.n_elem == 2 && "Check the model! It seems to be not compatible the linear model");

  return model(0) * x + model(1);
}

double pass::regression::predict_poly(const double &x, const arma::rowvec &model)
{
  assert(model.n_elem == 4 && "Check the model! It seems to be not compatible the polynomial model");

  return model(0) + model(1) * x + model(2) * std::pow(x, 2) + model(3) * std::pow(x, 3);
}
