#include "pass_bits/analyser/problem_evaluation_time.hpp"
#include <chrono> // assert

double pass::problem_evaluation_time(const pass::problem &problem)
{
  arma::vec particle(problem.dimension(), arma::fill::randu);
  arma::vec times(100);
  std::chrono::high_resolution_clock::time_point start, end;

  int runs = 0;

  while (runs < 200)
  {
    start = std::chrono::high_resolution_clock::now();
    problem.evaluate_normalised(particle);
    end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    times[runs % 100] = duration;
    ++runs;
  }

  return arma::median(times);
}
