#include "pass_bits/analyser/problem_evaluation_time.hpp"
#include <chrono> // assert

double pass::problem_evaluation_time(const pass::problem &problem)
{
  arma::vec particle(problem.dimension(), arma::fill::randu);
  arma::vec times(1000);
  std::chrono::high_resolution_clock::time_point start, end;

  // warm up
  for (arma::uword i = 0; i < 100000; ++i)
  {
    start = std::chrono::high_resolution_clock::now();
    problem.evaluate_normalised(particle);
    end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    times[i % 1000] = duration;
  }

  // do the evaluations
  for (arma::uword i = 0; i < 20000; ++i)
  {
    start = std::chrono::high_resolution_clock::now();
    problem.evaluate_normalised(particle);
    end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    times[i % 1000] = duration;
  }

  return arma::mean(times);
}
