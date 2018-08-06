#include "pass_bits/pagmo2_bindings.hpp"
#include <stdexcept>

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
    return pagmo::vector_double{wrapped_problem->evaluate_normalised({agent})};
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
    : optimiser(name),
      population_size(0) {}

pass::optimise_result pass::pagmo2::algorithm_adapter::optimise(
    const pass::problem &problem)
{
    pass::stopwatch stopwatch;
    stopwatch.start();

    pass::pagmo2::problem_adapter adapter{problem};
    pagmo::problem pagmo_problem{adapter};
    pagmo::population population{pagmo_problem, population_size};
    pagmo::algorithm algorithm = get_algorithm(problem);

    pass::optimise_result result{problem, acceptable_fitness_value};

    do
    {
        population = algorithm.evolve(population);
        result.normalised_agent = population.champion_x();
        result.fitness_value = population.champion_f()[0];
        result.duration = stopwatch.get_elapsed();
        result.iterations++;
        result.evaluations = pagmo_problem.get_fevals();
    } while (result.duration < maximal_duration &&
             result.iterations < maximal_iterations && !result.solved());

    return result;
}

// ----------------------------------------------
// CMAES
// ----------------------------------------------

pass::pagmo2::cmaes::cmaes() noexcept
    : algorithm_adapter("CMAES") {}

pagmo::algorithm pass::pagmo2::cmaes::get_algorithm(const pass::problem &) const
{
    return pagmo::algorithm{pagmo::cmaes{
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
        // Use -infinity because we have our own termination criteria check
        -std::numeric_limits<double>::infinity(),
        // xtol: stopping criteria on the f tolerance (default is 1e-6)
        // Use -infinity because we have our own termination criteria check
        -std::numeric_limits<double>::infinity(),
        // memory: when true the adapted parameters are not reset between successive
        // calls to the evolve method (default value: false)
        true,
        // force_bounds: when true the box bounds are enforced. The fitness will
        // never be called outside the bounds but the covariance matrix adaptation
        // mechanism will worsen (default value: false)
        true}};
}

// ----------------------------------------------
// differential evolution
// ----------------------------------------------

pass::pagmo2::differential_evolution::differential_evolution() noexcept
    : algorithm_adapter("differential_evolution") {}

pagmo::algorithm pass::pagmo2::differential_evolution::get_algorithm(const pass::problem &) const
{
    return pagmo::algorithm{pagmo::de{
        // gen: number of generations. (default value: 1)
        1,
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

// ----------------------------------------------
// simple genetic algorithm
// ----------------------------------------------

pass::pagmo2::simple_genetic_algorithm::simple_genetic_algorithm() noexcept
    : algorithm_adapter("simple_genetic_algorithm") {}

pagmo::algorithm pass::pagmo2::simple_genetic_algorithm::get_algorithm(const pass::problem &) const
{
    return pagmo::algorithm{pagmo::sga{}};
}

// ----------------------------------------------
// artifical bee colony
// ----------------------------------------------

pass::pagmo2::artifical_bee_colony::artifical_bee_colony() noexcept
    : algorithm_adapter("artifical_bee_colony") {}

pagmo::algorithm pass::pagmo2::artifical_bee_colony::get_algorithm(const pass::problem &) const
{
    return pagmo::algorithm{pagmo::bee_colony{}};
}

// ----------------------------------------------
// compass search
// ----------------------------------------------

const double pass::pagmo2::compass_search::initial_step_size = 0.3;
const double pass::pagmo2::compass_search::step_size_tolerance = 0.01;

pass::pagmo2::compass_search::compass_search() noexcept
    : algorithm_adapter("compass_search") {}

pagmo::algorithm pass::pagmo2::compass_search::get_algorithm(const pass::problem &problem) const
{
    return pagmo::algorithm{pagmo::compass_search{
        // max_fevals: maximum number of fitness evaluations
        static_cast<unsigned int>(maximal_iterations * problem.dimension() * 2),
        // start_range: start range (default is 0.1)
        initial_step_size,
        // stop_range: stop range (default is 0.01)
        step_size_tolerance}};
}
}
