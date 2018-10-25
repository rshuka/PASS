#include "pass_bits/analyser/openmp.hpp"
#include "pass_bits/problem/space_mission/gtoc1.hpp"
#include "pass_bits/helper/evaluation_time_stall.hpp"
#include "pass_bits/optimiser/parallel_swarm_search.hpp"
#include "pass_bits/optimiser/particle_swarm_optimisation.hpp"
#include "pass_bits/analyser/problem_evaluation_time.hpp"
#include "pass_bits/helper/regression.hpp"

bool pass::enable_openmp(const pass::problem &problem)
{
  // set random seed
  arma::arma_rng::set_seed_random();

  // define the maximum of runs
  arma::uword alg_runs = 2;

  // Array including all alg runtime, we want to test
  std::array<int, 30> repetitions = {1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15, 17, 20, 25, 28, 30, 33, 36,
                                     40, 45, 50, 60, 70, 80, 100, 120, 140, 160};

  // Output information
  std::cout << " =========================== Start openMP Analyse ========================= " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << " Your Problem:                " << problem.name << std::endl;
  std::cout << " Dimension:                   " << problem.dimension() << std::endl;

  std::cout << " ============================= Start Evaluation =========================== " << std::endl;

  double your_time = pass::problem_evaluation_time(problem);

  std::cout << "                                                                            " << std::endl;
  std::cout << " Evaluation time: " << your_time * 1e-6 << " microseconds." << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << " ============================= End Evaluation  ============================ " << std::endl;
  std::cout << "                                                                            " << std::endl;

  // Output information
  std::cout << " ============================= Start Trainining =========================== " << std::endl;
  std::cout << "                                                                            " << std::endl;

  arma::vec serial(alg_runs);
  arma::vec parallel(alg_runs);

  arma::mat summary(2, repetitions.size());

  pass::particle_swarm_optimisation algorithm_serial;
  algorithm_serial.maximal_duration = std::chrono::seconds(5);

  pass::parallel_swarm_search algorithm_parallel;
  algorithm_parallel.maximal_duration = std::chrono::seconds(5);

  std::srand(time(0));
  int count = 0;

  for (auto repetition : repetitions)
  {
    // Problem initialisation
    pass::gtoc1 test_problem;
    pass::evaluation_time_stall simulated_problem(test_problem);
    simulated_problem.repetitions = repetition;

    double ev_time = pass::problem_evaluation_time(simulated_problem);
    summary(0, count) = ev_time;

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

    summary(1, count) = arma::median(parallel) / arma::median(serial);
    count++;

    // load bar
    double temp_count = 30.0 / (count + 1);
    std::cout << " \r " << 100.0 / temp_count << " % completed." << std::flush;
  }

  std::cout << std::endl
            << std::endl
            << " Training completed successfully.\n"
            << std::flush;
  std::cout << "                                                                            " << std::endl;
  std::cout << " ===========================  End Training  =============================== " << std::endl;

  // Generating a regression object
  pass::regression r;

  // Getting the data for the model
  arma::rowvec x_values = summary.row(0);
  arma::rowvec y_values = summary.row(1);

  // Output information
  std::cout << " =========================== Start Building Models ======================== " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << " Building linear model for the training data.                               " << std::endl;
  std::cout << "                                                                            " << std::endl;

  // Generating the linear model

  arma::rowvec linear_model = r.linear_model(x_values, y_values);

  if (linear_model(2) >= 0.9)
  {
    std::cout << " Linear Model is suitable.                                                  " << std::endl;
    double predict_linear = r.predict_linear(your_time, linear_model);

    if (predict_linear > 0.9 * pass::number_of_threads())
    {
      predict_linear = 0.9 * pass::number_of_threads();
    }
    std::cout << "                                                                            " << std::endl;
    std::cout << " Your speedUp would be approximately: " << predict_linear << std::endl;

    if (predict_linear > pass::number_of_threads() / 2) // efficienty is more than 50 %
    {
      std::cout << "                                                                            " << std::endl;
      std::cout << " You should activate openMP!                                                " << std::endl;
      return true;
    }
    if (predict_linear < 1) // is efficienty is more than 50 %
    {
      std::cout << "                                                                            " << std::endl;
      std::cout << " You should NOT activate openMP!                                            " << std::endl;
      return false;
    }
    if (predict_linear > 1 && predict_linear < pass::number_of_threads() / 2) // efficienty is less than 50 %
    {
      std::cout << "                                                                            " << std::endl;
      std::cout << " You should decide yourself if to activate openMP or not                    " << std::endl;
    }
  }
  std::cout << " Linear Model is NOT suitable.                                              " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << " Finished building linear model                                            " << std::endl;

  std::cout << "                                                                            " << std::endl;
  std::cout << " Building polynomial model for the training data.                           " << std::endl;
  std::cout << "                                                                            " << std::endl;

  arma::rowvec poly_model = r.poly_model(x_values, y_values);

  double predict_poly = r.predict_poly(your_time, poly_model);

  if (predict_poly > 0.9 * pass::number_of_threads())
  {
    predict_poly = 0.9 * pass::number_of_threads();
  }
  std::cout << "                                                                            " << std::endl;
  std::cout << " Your speedUp would be approximately: " << predict_poly << std::endl;

  if (predict_poly > pass::number_of_threads() / 2) // is efficienty is more than 50 %
  {
    std::cout << "                                                                            " << std::endl;
    std::cout << " You should activate openMP!                                                " << std::endl;
    return true;
  }
  if (predict_poly < 1) // is efficienty is more than 50 %
  {
    std::cout << "                                                                            " << std::endl;
    std::cout << " You should NOT activate openMP!                                            " << std::endl;
    return false;
  }
  if (predict_poly > 1 && predict_poly < pass::number_of_threads() / 2) // is efficienty is more than 50 %
  {
    std::cout << "                                                                            " << std::endl;
    std::cout << " You should decide yourself if to activate openMP or not                    " << std::endl;
  }

  std::cout << " ========================= Done Building Models  ========================== " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << " =========================  Done openMP Analyse  ========================== " << std::endl;

  return true;
}
