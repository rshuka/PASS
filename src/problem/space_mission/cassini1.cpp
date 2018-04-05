#include "pass_bits/problem/space_mission/cassini1.hpp"
#include "pass_bits/helper/astro_problems/mga.hpp"

pass::cassini1::cassini1()
    : problem({-1000, 30, 100, 30, 400, 1000},
              {0, 400, 470, 400, 2000, 6000},
              "Cassini1") {}

double pass::cassini1::evaluate(const arma::vec &agent) const
{
  assert(agent.n_elem == dimension() &&
         "`agent` has incompatible dimension");

  std::vector<double> x(dimension());
  for (arma::uword n = 0; n < agent.n_elem; n++)
  {
    x[n] = agent[n];
  }
  std::vector<double> rp(dimension());

  const int CASSINI_DIM = 6;
  std::vector<double> Delta_V(CASSINI_DIM);
  rp.resize(CASSINI_DIM - 2);
  std::vector<double> t(CASSINI_DIM);
  mgaproblem problem;

  //Filling up the problem parameters
  problem.type = 1; // type 1 is Cassini 1; compare with mga_dsm.cpp

  int sequence_[CASSINI_DIM] = {3, 2, 2, 3, 5, 6}; // sequence of planets
  std::vector<int> sequence(CASSINI_DIM);
  problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_ + CASSINI_DIM);

  const int rev_[CASSINI_DIM] = {0, 0, 0, 0, 0, 0}; // sequence of clockwise legs
  std::vector<int> rev(CASSINI_DIM);
  problem.rev_flag.insert(problem.rev_flag.begin(), rev_, rev_ + CASSINI_DIM);

  problem.e = 0.98;     // Final orbit eccentricity
  problem.rp = 108950;  // Final orbit pericenter
  problem.DVlaunch = 0; // Launcher DV

  double obj = 0;

  MGA(x, problem, rp, Delta_V, obj);

  return obj;
}
