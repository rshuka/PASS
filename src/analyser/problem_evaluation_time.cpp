#include "pass_bits/analyser/problem_evaluation_time.hpp"
#include <chrono>

double pass::problem_evaluation_time(const pass::problem &problem)
{
  arma::vec times(1000);
  std::chrono::high_resolution_clock::time_point start, end;

  // warm up
  for (arma::uword i = 0; i < 10000; ++i)
  {
    arma::vec particle(problem.dimension(), arma::fill::randu);
    assert(particle.n_elem == problem.dimension() &&
           "`particle` has incompatible dimension");

    start = std::chrono::high_resolution_clock::now();
    problem.evaluate_normalised(particle);
    end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    times[i % 100] = duration;
  }

  // do the evaluations
  for (arma::uword i = 0; i < 2000; ++i)
  {
    arma::vec particle(problem.dimension(), arma::fill::randu);
    assert(particle.n_elem == problem.dimension() &&
           "`particle` has incompatible dimension");

    start = std::chrono::high_resolution_clock::now();
    problem.evaluate_normalised(particle);
    end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    times[i % 1000] = duration;
  }

  return arma::mean(times);
}
