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
  arma::uword alg_runs = 3;

  // Array including all alg runtime, we want to test
  std::array<int, 30> repetitions = {1, 2, 3, 6, 8, 10, 14, 17, 20, 25, 28, 30, 33, 36, 40,
                                     45, 50, 60, 70, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 300};

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
    pass::gtoc1 problem;
    pass::evaluation_time_stall simulated_problem(problem);
    simulated_problem.repetitions = repetition;

    double t = pass::problem_evaluation_time(simulated_problem);
    std::cout << "Time: " << t * 1e-6 << std::endl;

    // Do the evaluation for serial and parallel for all the evaluations values
    for (arma::uword serial_run = 0; serial_run < alg_runs; ++serial_run)
    {
      auto serial_alg = algorithm_serial.optimise(simulated_problem);
      serial(serial_run) = serial_alg.evaluations;
    }

    summary(count, 0) = arma::median(serial);

    for (arma::uword parallel_run = 0; parallel_run < alg_runs; ++parallel_run)
    {
      auto parallel_alg = algorithm_parallel.optimise(simulated_problem);
      parallel(parallel_run) = parallel_alg.evaluations;
    }

    summary(count, 1) = arma::median(parallel);

    std::cout << "Speedup: " << (arma::median(parallel) / arma::median(serial)) << std::endl;
    std::cout << "---------------------------" << std::endl;
    count++;
  }

  return true;
}
