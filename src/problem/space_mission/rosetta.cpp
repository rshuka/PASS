#include "pass_bits/problem/space_mission/rosetta.hpp"
#include "pass_bits/helper/astro_problems/mga_dsm.hpp"

pass::rosetta::rosetta()
    : problem({1460, 3, 0, 0, 300, 150, 150, 300, 700, 0.01, 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.05, 1.05, -arma::datum::pi, -arma::datum::pi, -arma::datum::pi, -arma::datum::pi},
              {1825, 5, 1, 1, 500, 800, 800, 800, 1850, 0.9, 0.9, 0.9, 0.9, 0.9, 9, 9, 9, 9, arma::datum::pi, arma::datum::pi, arma::datum::pi, arma::datum::pi},
              "Rosetta") {}

double pass::rosetta::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");

  std::vector<double> x(dimension());
  for (arma::uword n = 0; n < agent.n_elem; n++)
  {
    x[n] = agent[n];
  }

  mgadsmproblem problem;

  int sequence_[6] = {3, 3, 4, 3, 3, 10}; // sequence of planets
  problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_ + 6);
  problem.type = 2; // type 2 means rosetta; compare with mga_dsm.cpp
  problem.asteroid.keplerian[0] = 3.50294972836275;
  problem.asteroid.keplerian[1] = 0.6319356;
  problem.asteroid.keplerian[2] = 7.12723;
  problem.asteroid.keplerian[3] = 50.92302;
  problem.asteroid.keplerian[4] = 11.36788;
  problem.asteroid.keplerian[5] = 0.0;
  problem.asteroid.epoch = 52504.23754000012;
  problem.asteroid.mu = 0.0;

  //Allocate temporary memory for MGA_DSM
  problem.r = std::vector<double *>(6);
  problem.v = std::vector<double *>(6);
  problem.DV = std::vector<double>(6 + 1);

  for (int i = 0; i < 6; i++)
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

  //Free temporary memory for MGA_DSM
  for (int i = 0; i < 6; i++)
  {
    delete[] problem.r[i];
    delete[] problem.v[i];
  }
  problem.r.clear();
  problem.v.clear();

  return obj;
}
