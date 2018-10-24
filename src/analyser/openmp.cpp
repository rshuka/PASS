#include "pass_bits/analyser/openmp.hpp"
#include "pass_bits/problem/optimisation_benchmark/styblinski_tang_function.hpp"
#include "pass_bits/helper/evaluation_time_stall.hpp"
#include "pass_bits/optimiser/parallel_swarm_search.hpp"
#include "pass_bits/optimiser/particle_swarm_optimisation.hpp"
#include <array>

bool pass::enable_openmp()
{
  // set random seed
  arma::arma_rng::set_seed_random();

  // define the maximum of runs
  arma::uword alg_runs = 5;

  // Array including all alg runtime, we want to test
  std::array<int, 15> repetitions = {1, 10, 50, 90, 140, 190, 200, 300, 900, 1200, 1800, 2000, 2400, 2800, 3000};

  arma::vec serial(alg_runs);
  arma::vec parallel(alg_runs);

  arma::mat summary(repetitions.size(), 2);

  pass::particle_swarm_optimisation algorithm_serial;
  algorithm_serial.maximal_duration = std::chrono::seconds(5);

  pass::parallel_swarm_search algorithm_parallel;
  algorithm_parallel.maximal_duration = std::chrono::seconds(5);

  int count = 0;

  for (auto repetition : repetitions)
  {
    std::cout << "Repetition: " << repetition << std::endl;

    // Problem initialisation
    pass::styblinski_tang_function problem(10);
    pass::evaluation_time_stall simulated_problem(problem);
    simulated_problem.repetitions = repetition;

    // Do the evaluation for serial and parallel for all the evaluations values

    for (arma::uword serial_run = 0; serial_run < alg_runs; ++serial_run)
    {
      auto serial_alg = algorithm_serial.optimise(simulated_problem);
      serial(serial_run) = serial_alg.iterations;
    }

    std::cout << "Serial list \n: " << serial << std::endl;

    summary(count, 0) = arma::median(serial);

    std::cout << "Summary \n: " << summary << std::endl;

    for (arma::uword parallel_run = 0; parallel_run < alg_runs; ++parallel_run)
    {
      auto parallel_alg = algorithm_parallel.optimise(simulated_problem);
      parallel(parallel_run) = parallel_alg.iterations;
    }
    std::cout << "Parallel list \n: " << parallel << std::endl;

    summary(count, 1) = arma::median(parallel);
    std::cout << "Summary \n: " << summary << std::endl;

    std::cout << "Speedup: " << (arma::median(parallel) / arma::median(serial)) << std::endl;
    std::cout << "---------------------------" << std::endl;
    count++;
  }

  return true;
}
