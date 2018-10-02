#pragma once
#include "pass_bits/problem.hpp"
#include "pass_bits/optimiser.hpp"

namespace pass
{
/**
 * Search the best parameters for the parallel swarm search algorithm
 * and saves the parameters in a file which can be than loaded
 *
 * Following parameter of the parallel swarm search are investigated
 * - Swarm Size
 * - Neighbourhood Probability
 * - Inertia
 * - Cognitive Acceleration = Social Acceleration
 *
 * @parameters
 * - problem: the given problem
 * - Benchmark? "True" or "False"
 */

void search_parameters(const pass::problem &problem, const bool benchmark);

/**
 * Evaluate the Problem
 *
 * The Problem will be evaluated 10 times
 * Set the global variable to change that
 * pass::parameter_setting_number_of_runs
 */

arma::mat parameter_evaluate(pass::optimiser &optimiser, const pass::problem &problems);

/**
 * Compare the two segments with each other
 *
 * @return
 * - 1 if the first segmet is the best
 * - 2 if the second segmet is the best
 */

arma::uword compare_segments(const arma::mat first_segment_runtimes, const arma::mat second_segment_runtimes);

} // namespace pass
