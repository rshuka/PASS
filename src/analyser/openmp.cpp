#include "pass_bits/analyser/openmp.hpp"
#include "pass_bits/problem/space_mission/gtoc1.hpp"
#include "pass_bits/helper/evaluation_time_stall.hpp"
#include "pass_bits/optimiser/parallel_swarm_search.hpp"
#include "pass_bits/optimiser/particle_swarm_optimisation.hpp"
#include "pass_bits/analyser/problem_evaluation_time.hpp"
#include <array>

bool pass::enable_openmp()
{
  // set random seed
  arma::arma_rng::set_seed_random();

  // define the maximum of runs
  arma::uword alg_runs = 5;

  // Array including all alg runtime, we want to test
  std::array<int, 15> repetitions = {1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140};

  arma::vec serial(alg_runs);
  arma::vec parallel(alg_runs);

  arma::mat summary(repetitions.size(), 2);

  pass::particle_swarm_optimisation algorithm_serial;
  //algorithm_serial.maximal_duration = std::chrono::seconds(10);
  algorithm_serial.maximal_iterations = 100;

  pass::parallel_swarm_search algorithm_parallel;
  //algorithm_parallel.maximal_duration = std::chrono::seconds(10);
  algorithm_parallel.maximal_iterations = 100;

  int count = 0;

  for (auto repetition : repetitions)
  {
    std::cout << "Repetition: " << repetition << std::endl;

    // Problem initialisation
    pass::gtoc1 problem;
    pass::evaluation_time_stall simulated_problem(problem);
    simulated_problem.repetitions = repetition;

    double t = pass::problem_evaluation_time(simulated_problem);

    std::cout << "Time: " << t * 1e-6 << std::endl;
    // Do the evaluation for serial and parallel for all the evaluations values

    for (arma::uword serial_run = 0; serial_run < alg_runs; ++serial_run)
    {
      auto serial_alg = algorithm_serial.optimise(simulated_problem);
      //serial(serial_run) = serial_alg.iterations;
      serial(serial_run) = serial_alg.duration.count();
    }

    std::cout << "Serial list \n: " << serial << std::endl;

    summary(count, 0) = arma::median(serial);

    std::cout << "Summary \n: " << summary << std::endl;

    for (arma::uword parallel_run = 0; parallel_run < alg_runs; ++parallel_run)
    {
      auto parallel_alg = algorithm_parallel.optimise(simulated_problem);
      //parallel(parallel_run) = parallel_alg.iterations;
      parallel(parallel_run) = parallel_alg.duration.count();
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
