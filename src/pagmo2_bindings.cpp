#include "pass_bits/pagmo2_bindings.hpp"
#include <algorithm>
#include <stdexcept>
#include <cmath>

// ----------------------------------------------
// problem adapter
// ----------------------------------------------

pass::pagmo2::problem_adapter::problem_adapter()
    : wrapped_problem(nullptr)
{
    throw std::domain_error("pass::pagmo2::problem_adapter is not really default constructible");
}

pass::pagmo2::problem_adapter::problem_adapter(const pass::problem &wrapped_problem)
    : wrapped_problem(&wrapped_problem)
{
}

pagmo::vector_double pass::pagmo2::problem_adapter::fitness(const pagmo::vector_double &agent) const
{
    arma::vec v = {agent};
    if (arma::any(v < 0) || arma::any(v > 1))
    {
        return {std::numeric_limits<double>::infinity()};
    }
    return pagmo::vector_double{wrapped_problem->evaluate_normalised(v)};
}

std::pair<pagmo::vector_double, pagmo::vector_double> pass::pagmo2::problem_adapter::get_bounds() const
{
    return {pagmo::vector_double(wrapped_problem->dimension(), 0),
            pagmo::vector_double(wrapped_problem->dimension(), 1)};
}

// ----------------------------------------------
// algorithm adapter
// ----------------------------------------------

pass::pagmo2::algorithm_adapter::algorithm_adapter(const std::string &name) noexcept
    : optimiser(name) {}

pass::optimise_result pass::pagmo2::algorithm_adapter::optimise(
    const pass::problem &problem)
{
    assert(maximal_iterations == std::numeric_limits<arma::uword>::max() &&
           "The pagmo optimisers don't respect the `maximal_iterations` termination criterion");
    assert(acceptable_fitness_value == -std::numeric_limits<double>::infinity() &&
           "The pagmo optimisers don't respect the `acceptable_fitness_value` termination criterion");
    assert(maximal_duration == std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::duration::max()) &&
           "The pagmo optimisers don't respect the `maximal_duration` termination criterion");

    pass::stopwatch stopwatch;
    stopwatch.start();

    pass::pagmo2::problem_adapter adapter{problem};
    pagmo::problem pagmo_problem{adapter};
    pagmo::population population{pagmo_problem, calculate_population_size(problem)};
    pagmo::algorithm algorithm = get_algorithm(problem);

    pass::optimise_result result{problem, acceptable_fitness_value};
    population = algorithm.evolve(population);
    result.normalised_agent = population.champion_x();
    result.fitness_value = population.champion_f()[0];
    result.duration = stopwatch.get_elapsed();
    result.iterations = calculate_iterations(problem);
    result.evaluations = population.get_problem().get_fevals();

    return result;
}

// ----------------------------------------------
// CMAES
// ----------------------------------------------

pass::pagmo2::cmaes::cmaes() noexcept
    : algorithm_adapter("CMAES") {}

pagmo::algorithm pass::pagmo2::cmaes::get_algorithm(const pass::problem &problem) const
{
    return pagmo::algorithm{pagmo::cmaes{
        // gen: number of generations (default value: 1)
        static_cast<unsigned>(calculate_iterations(problem)),
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
        // Use -infinity because we have our own termination criteria check
        -std::numeric_limits<double>::infinity(),
        // xtol: stopping criteria on the f tolerance (default is 1e-6)
        // Use -infinity because we have our own termination criteria check
        -std::numeric_limits<double>::infinity()}};
}

arma::uword pass::pagmo2::cmaes::calculate_iterations(const pass::problem &problem) const
{
    return minimal_evaluations / calculate_population_size(problem) +
           (minimal_evaluations % calculate_population_size(problem) ? 1 : 0);
}

arma::uword pass::pagmo2::cmaes::calculate_population_size(const pass::problem &problem) const
{
    return 4 + 3 * static_cast<arma::uword>(log(problem.dimension()));
}

// ----------------------------------------------
// differential evolution
// ----------------------------------------------

pass::pagmo2::differential_evolution::differential_evolution() noexcept
    : algorithm_adapter("differential_evolution") {}

pagmo::algorithm pass::pagmo2::differential_evolution::get_algorithm(const pass::problem &problem) const
{
    return pagmo::algorithm{pagmo::de{
        // gen: number of generations. (default value: 1)
        static_cast<unsigned>(calculate_iterations(problem)),
        // F: weight coefficient (dafault value is 0.8)
        0.8,
        // CR: crossover probability (dafault value is 0.9)
        0.9,
        // variant: mutation variant (dafault variant is 2: /rand/1/exp)
        2,
        // ftol: stopping criteria on the x tolerance (default is 1e-6)
        // Use -infinity because we have our own termination criteria check
        -std::numeric_limits<double>::infinity(),
        // xtol: stopping criteria on the f tolerance (default is 1e-6)
        // Use -infinity because we have our own termination criteria check
        -std::numeric_limits<double>::infinity()}};
}

arma::uword pass::pagmo2::differential_evolution::calculate_iterations(const pass::problem &problem) const
{
    return minimal_evaluations / calculate_population_size(problem) +
           (minimal_evaluations % calculate_population_size(problem) ? 1 : 0);
}

arma::uword pass::pagmo2::differential_evolution::calculate_population_size(const pass::problem &) const
{
    return 5;
}

// ----------------------------------------------
// simple genetic algorithm
// ----------------------------------------------

pass::pagmo2::simple_genetic_algorithm::simple_genetic_algorithm() noexcept
    : algorithm_adapter("simple_genetic_algorithm") {}

pagmo::algorithm pass::pagmo2::simple_genetic_algorithm::get_algorithm(const pass::problem &problem) const
{
    return pagmo::algorithm{pagmo::sga{
        // gen:
        static_cast<unsigned>(calculate_iterations(problem))}};
}

arma::uword pass::pagmo2::simple_genetic_algorithm::calculate_iterations(const pass::problem &problem) const
{
    return minimal_evaluations / calculate_population_size(problem) +
           (minimal_evaluations % calculate_population_size(problem) ? 1 : 0);
}

arma::uword pass::pagmo2::simple_genetic_algorithm::calculate_population_size(const pass::problem &) const
{
    return 5;
}

// ----------------------------------------------
// artifical bee colony
// ----------------------------------------------

pass::pagmo2::artifical_bee_colony::artifical_bee_colony() noexcept
    : algorithm_adapter("artifical_bee_colony") {}

pagmo::algorithm pass::pagmo2::artifical_bee_colony::get_algorithm(const pass::problem &problem) const
{
    return pagmo::algorithm{pagmo::bee_colony{
        // gen:
        static_cast<unsigned>(calculate_iterations(problem))}};
}

arma::uword pass::pagmo2::artifical_bee_colony::calculate_iterations(const pass::problem &problem) const
{
    arma::uword evaluations_per_iteration = 2 * calculate_population_size(problem);
    return minimal_evaluations / evaluations_per_iteration +
           (minimal_evaluations % evaluations_per_iteration ? 1 : 0);
}

arma::uword pass::pagmo2::artifical_bee_colony::calculate_population_size(const pass::problem &) const
{
    return 5;
}

// ----------------------------------------------
// compass search
// ----------------------------------------------

const double pass::pagmo2::compass_search::initial_step_size = 0.3;
const double pass::pagmo2::compass_search::step_size_tolerance = 0.01;

pass::pagmo2::compass_search::compass_search() noexcept
    : algorithm_adapter("compass_search") {}

pagmo::algorithm pass::pagmo2::compass_search::get_algorithm(const pass::problem &) const
{
    return pagmo::algorithm{pagmo::compass_search{
        // max_fevals: maximum number of fitness evaluations
        static_cast<unsigned>(minimal_evaluations),
        // start_range: start range (default is 0.1)
        initial_step_size,
        // stop_range: stop range (default is 0.01)
        0}};
}

arma::uword pass::pagmo2::compass_search::calculate_iterations(const pass::problem &problem) const
{
    return minimal_evaluations / (2 * problem.dimension());
}

arma::uword pass::pagmo2::compass_search::calculate_population_size(const pass::problem &) const
{
    return 1;
}

// ----------------------------------------------
// simulated annealing
// ----------------------------------------------

pass::pagmo2::simulated_annealing::simulated_annealing() noexcept
    : algorithm_adapter("simulated_annealing") {}

pagmo::algorithm pass::pagmo2::simulated_annealing::get_algorithm(const pass::problem &problem) const
{
    return pagmo::algorithm{pagmo::simulated_annealing{
        // Ts: starting temperature (default is 10.0)
        10.0,
        // Tf: final temperature (default is 0.1)
        0.1,
        // n_T_adj: number of temperature adjustments in the annealing schedule
        // (default is 10)
        static_cast<unsigned>(calculate_iterations(problem)),
        // n_range_adj: number of adjustments of the search range performed at a
        // constant temperature (default is 1)
        std::max<unsigned>(100, 5 * problem.dimension()),
        // bin_size: number of mutations that are used to compute the acceptance
        // rate (default is 20)
        20,
        // start_range: starting range for mutating the decision vector
        // (default is 1)
        1.0}};
}

arma::uword pass::pagmo2::simulated_annealing::calculate_iterations(const pass::problem &problem) const
{
    arma::uword evaluations_per_iteration = 20 * std::max<arma::uword>(100, 5 * problem.dimension()) * problem.dimension();
    return minimal_evaluations / evaluations_per_iteration +
           (minimal_evaluations % evaluations_per_iteration ? 1 : 0);
}

arma::uword pass::pagmo2::simulated_annealing::calculate_population_size(const pass::problem &) const
{
    return 1;
}
