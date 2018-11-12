#pragma once

// Armadillo
#define ARMA_DONT_PRINT_ERRORS
#define ARMA_USE_CXX11
#include <armadillo>

// MPI support must be added via CMake, to ensure that we also link against it.
// Therefore, CMake will decide whether SUPPORT_MPI is to be defined or not.
#cmakedefine SUPPORT_MPI

// OpenMP support must be added via CMake, to ensure that we also link against it.
// Therefore, CMake will decide whether SUPPORT_OPENMP is to be defined or not.
#cmakedefine SUPPORT_OPENMP

// The maximal number of threads to be supported by PASS.
// Larger values may result in a greater start up time and decrese efficency.
// In case `MAXIMAL_NUMBER_OF_THREADS` was not defined before, we fall back to the value below, determined via CMake.
#if !defined(MAXIMAL_NUMBER_OF_THREADS)

#define MAXIMAL_NUMBER_OF_THREADS @MAXIMAL_NUMBER_OF_THREADS@

#endif

namespace pass
{
  /**
  * Global variables used for evaluations
  * @is_verbose: analyse the behaviour of the algorithm
  * Is initialized to `false`.
  */
  extern bool is_verbose;

  /**
  * Global variables used for evaluations
  * @global_number_of_runs: to help saving raw data
  * Starts from '1'
  */
  extern int global_number_of_runs;
  /**
  * Global variables used for evaluations
  * @precision: the precision to stop algorithms and set restart
  * Is initialized to `1e-06`.
  */
  extern double precision;
  /**
  * Global variables used for evaluations
  * @parameter_setting_number_of_runs: number of the runs
  * Is initialized to `10`.
  */
  extern arma::uword parameter_setting_number_of_runs;

  int thread_number();
  int number_of_threads();
  int node_rank();
  int number_of_nodes();
}
