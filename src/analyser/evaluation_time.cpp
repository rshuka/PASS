#include "pass_bits/analyser/problem_time.hpp"
#include <chrono> // assert
#include <thread>

double pass::evaluation_time(const pass::problem &problem)
{
  arma::vec particle(problem.dimension(), arma::fill::randu);
  arma::vec times(1000);
  std::chrono::high_resolution_clock::time_point t1, t2;

  int runs = 0;
  int count = 0;

  while (runs < 5000)
  {
    t1 = std::chrono::high_resolution_clock::now();
    problem.evaluate_normalised(particle);
    t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();

    times[count] = duration;

    count = (count + 1) % 1000;
    runs++;
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }

  return arma::mean(times);
}
