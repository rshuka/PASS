#include "pass_bits/optimiser/cmaes.hpp"
#include "pass_bits/helper/pagmo2/problem_adapter.hpp"
#include <pagmo/pagmo.hpp>

pass::cmaes::cmaes() noexcept
    : optimiser("CMAES"),
      population_size(0) {}

pass::optimise_result pass::cmaes::optimise(
    const pass::problem &problem)
{
  assert(population_size >= 5 && "");
  pass::stopwatch stopwatch;
  stopwatch.start();

  pass::pagmo2::problem_adapter adapter{problem};
  pagmo::problem pagmo_problem{adapter};
  pagmo::population population{pagmo_problem, population_size};
  pagmo::algorithm cmaes{pagmo::cmaes{
    // gen: number of generations (default value: 1)
    1,
    // cc: backward time horizon for the evolution path (default value: -1)
    -1,
    // cs: makes partly up for the small variance loss in case the indicator is
    // zero (default value: -1)
    -1,
    // c1: learning rate for the rank-one update of the covariance matrix
    // (default value: -1)
    -1,
    // cmu: learning rate for the rank- μ update of the covariance matrix
    // (default value: -1)
    -1,
    // sigma0: initial step-size (default value: 0.5)
    0.5,
    // ftol: stopping criteria on the x tolerance (default is 1e-6)
    1e-6,
    // xtol: stopping criteria on the f tolerance (default is 1e-6)
    1e-6,
    // memory: when true the adapted parameters are not reset between successive
    // calls to the evolve method (default value: false)
    true,
    // force_bounds: when true the box bounds are enforced. The fitness will
    // never be called outside the bounds but the covariance matrix adaptation
    // mechanism will worsen (default value: false)
    true
  }};

  pass::optimise_result result{problem, acceptable_fitness_value};

  do {
    population = cmaes.evolve(population);
    result.normalised_agent = population.champion_x();
    result.best_found_value = population.champion_f()[0];
    result.duration = stopwatch.get_elapsed();
    result.iterations++;
    result.evaluations = pagmo_problem.get_fevals();
  }   while (result.duration < maximal_duration &&
         result.iterations < maximal_iterations && !result.solved());

  return result;
}
