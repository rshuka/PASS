#include "pass_bits/analyser/evaluation_time.hpp"
#include <chrono> // assert
#include <thread> // sleep_for

double pass::evaluation_time(const pass::problem &problem)
{
  arma::vec particle(problem.dimension(), arma::fill::randu);
  arma::vec times(1000);
  std::chrono::high_resolution_clock::time_point start, end;

  int runs = 0;

  while (runs < 5000)
  {
    start = std::chrono::high_resolution_clock::now();
    problem.evaluate_normalised(particle);
    end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    times[runs % 1000] = duration;
    ++runs;

    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }

  return arma::median(times);
}
