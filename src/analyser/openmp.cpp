#include "pass_bits/analyser/openmp.hpp"
#include "pass_bits/problem/space_mission/gtoc1.hpp"
#include "pass_bits/problem/optimisation_benchmark/ackley_function.hpp"
#include "pass_bits/helper/evaluation_time_stall.hpp"
#include "pass_bits/optimiser/parallel_swarm_search.hpp"
#include "pass_bits/optimiser/particle_swarm_optimisation.hpp"
#include "pass_bits/analyser/problem_evaluation_time.hpp"
#include "pass_bits/helper/random.hpp"

bool pass::enable_openmp(const pass::problem &problem, const int &training)
{
  // set random seed
  arma::arma_rng::set_seed_random();

  // min and max of training data
  int min = 10;
  int max = 2000;

  // define the maximum of runs
  arma::uword alg_runs = 2;

  arma::rowvec data_normalised(training, arma::fill::randu);

  arma::rowvec data(training);

  for (int i = 0; i < training; i++)
  {
    data(i) = data_normalised(i) * max + min;
  }

  arma::vec serial(alg_runs);
  arma::vec parallel(alg_runs);

  arma::mat summary(data.size(), 2);

  pass::particle_swarm_optimisation algorithm_serial;
  algorithm_serial.maximal_duration = std::chrono::seconds(5);

  pass::parallel_swarm_search algorithm_parallel;
  algorithm_parallel.maximal_duration = std::chrono::seconds(5);

  int count = 0;

  for (auto repetition : data)
  {
    // Problem initialisation
    pass::ackley_function problem(50);
    pass::evaluation_time_stall simulated_problem(problem);
    simulated_problem.repetitions = repetition;

    double ev_time = pass::problem_evaluation_time(simulated_problem);
    summary(count, 0) = ev_time;

    // Do the evaluation for serial and parallel for all the evaluations values
    for (arma::uword serial_run = 0; serial_run < alg_runs; ++serial_run)
    {
      auto serial_alg = algorithm_serial.optimise(simulated_problem);
      serial(serial_run) = serial_alg.evaluations;
    }

    for (arma::uword parallel_run = 0; parallel_run < alg_runs; ++parallel_run)
    {
      auto parallel_alg = algorithm_parallel.optimise(simulated_problem);
      parallel(parallel_run) = parallel_alg.evaluations;
    }

    summary(count, 1) = arma::median(parallel) / arma::median(serial);

    count++;
  }

  std::cout << "Summary: " << summary << std::endl;

  return true;
}
