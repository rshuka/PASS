#include "pass_bits/config.hpp"

// C++ standard library and MPI
#if defined(SUPPORT_MPI)
#include <mpi.h>
#include <stdexcept>
#endif

#if defined(SUPPORT_OPENMP)
#include <omp.h>
#endif

/**
 * If parallelization is activated some useful variables are provided
 */
namespace pass
{
/**
  * Global variables used for evaluations
  * @is_verbose: analyse the behaviour of the algorithm
  *
  * @number_of_runs: to help saving raw data
  *
  * @precision: the precision to stop algorithms and set restart
  */
bool is_verbose(false);
int global_number_of_runs(1);
double precision(1e-06);
arma::uword parameter_setting_number_of_runs(10);

/**
 * Use OpenMP
 * @number_of_threads: returns the number of the threads
 * @thread_number: returns which number the thread have
 */
int number_of_threads()
{
#if defined(SUPPORT_OPENMP)
  return omp_get_max_threads();
#else
  return 1;
#endif
}

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
 * @node_rank: returns which number a.k.a rank the node have
 */
int number_of_nodes()
{
#if defined(SUPPORT_MPI)
  int numberOfNodes;
  MPI_Comm_size(MPI_COMM_WORLD, &numberOfNodes);

  if (numberOfNodes < 0)
  {
    throw std::runtime_error(
        "number_of_nodes: Please check your MPI installation, as we got a "
        "negative number of nodes.");
  }

  return numberOfNodes;
#else
  return 1;
#endif
}

int node_rank()
{
#if defined(SUPPORT_MPI)
  int nodeRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeRank);

  if (nodeRank < 0)
  {
    throw std::runtime_error(
        "node_rank: Please check your MPI installation, as we got a negative "
        "node rank.");
  }

  return nodeRank;
#else
  return 0;
#endif
}

} // namespace pass
