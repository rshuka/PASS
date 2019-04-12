#include "pass_bits/config.hpp"

// C++ standard library and MPI
#if defined(SUPPORT_MPI)
#include <mpi.h>
#include <stdexcept>
#endif

#if defined(SUPPORT_OPENMP)
#include <omp.h>
#endif

namespace pass
{

/**
  * Global variables used for evaluations
  * @is_verbose: analyse the behaviour of the algorithm
  * Is initialized to `false`.
  */
bool is_verbose(false);

/**
  * Global variables used for evaluations
  * @global_number_of_runs: to help saving raw data
  * Starts from '1'
  */
int global_number_of_runs(1);

/**
  * Global variables used for evaluations
  * @precision: the precision to stop algorithms and set restart
  * Is initialized to `1e-06`.
  */
double precision(1e-06);

/**
  * Global variables used for evaluations
  * @parameter_setting_number_of_runs: number of the runs
  * Is initialized to `10`.
  */
arma::uword parameter_setting_number_of_runs(30);

/**
 * Use OpenMP
 * @number_of_threads: returns the number of the threads
 */
int number_of_threads()
{
#if defined(SUPPORT_OPENMP)
  return omp_get_max_threads();
#else
  return 1;
#endif
}

/**
 * Use OpenMP
 * @thread_number: returns which number the thread have
 */
int thread_number()
{
#if defined(SUPPORT_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}

/**
 * Use MPI
 * @number_of_nodes: returns the number of the nodes
 */
int number_of_nodes()
{
#if defined(SUPPORT_MPI)
  int number_of_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes);

  if (number_of_nodes < 0)
  {
    throw std::runtime_error(
        "number_of_nodes: Please check your MPI installation, as we got a "
        "negative number of nodes.");
  }

  return number_of_nodes;
#else
  return 1;
#endif
}

/**
 * Use MPI
 * @node_rank: returns which number a.k.a rank the node have
 */
int node_rank()
{
#if defined(SUPPORT_MPI)
  int node_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &node_rank);

  if (node_rank < 0)
  {
    throw std::runtime_error(
        "node_rank: Please check your MPI installation, as we got a negative "
        "node rank.");
  }

  return node_rank;
#else
  return 0;
#endif
}

} // namespace pass
