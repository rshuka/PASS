#include "pass_bits/problem/messenger_full.hpp"
#include "pass_bits/helper/gtoc1/mga_dsm.hpp"

pass::messenger_full::messenger_full()
    : problem({1900, 2.5, 0, 0, 100, 100, 100, 100, 100, 100, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.05, 1.05, 1.05, -arma::datum::pi, -arma::datum::pi, -arma::datum::pi, -arma::datum::pi, -arma::datum::pi},
              {2300, 4.05, 1, 1, 500, 500, 500, 500, 500, 600, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 6, 6, 6, 6, 6, arma::datum::pi, arma::datum::pi, arma::datum::pi, arma::datum::pi, arma::datum::pi},
              "Messenger_Full") {}

double pass::messenger_full::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");

  std::vector<double> x(dimension());
  for (arma::uword n = 0; n < agent.n_elem; n++)
  {
    x[n] = agent[n];
  }

  mgadsmproblem problem;

  int sequence_[7] = {3, 2, 2, 1, 1, 1, 1};
  // sequence of planets
  problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_ + 7);
  const int orbit_insertion = 0; // From: mga.h
  problem.type = orbit_insertion;
  problem.e = 0.704;
  problem.rp = 2640.0;

  //Memory allocation
  problem.r = std::vector<double *>(7);
  problem.v = std::vector<double *>(7);
  problem.DV = std::vector<double>(7 + 1);
  for (int i = 0; i < 7; i++)
  {
    problem.r[i] = new double[3];
    problem.v[i] = new double[3];
  }

  double obj = 0;

  MGA_DSM(
      /* INPUT values: */
      x,
      problem,

      /* OUTPUT values: */
      obj);

  //Memory release
  for (int i = 0; i < 7; i++)
  {
    delete[] problem.r[i];
    delete[] problem.v[i];
  }
  problem.r.clear();
  problem.v.clear();

  return obj;
}
