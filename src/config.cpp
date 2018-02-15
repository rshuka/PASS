#include "pass_bits/config.hpp"

// C++ standard library and MPI
#if defined(SUPPORT_MPI)
#include <mpi.h>
#include <stdexcept>
#endif

/**
 * If parallelization is activated some useful variables are provided
 */
namespace pass {
/**
 * Use OpenMP
 * @numberOfThreads: returns the number of the threads
 * @threadNumber: returns which number the thread have
 */
arma::uword numberOfThreads() {
#if defined(SUPPORT_OPENMP)
  return omp_get_max_threads();
#else
  return 1;
#endif
}

arma::uword threadNumber() {
#if defined(SUPPORT_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}

/**
 * Use MPI
 * @numberOfNodes: returns the number of the nodes
 * @nodeRank: returns which number a.k.a rank the node have
 */
arma::uword numberOfNodes() {
#if defined(SUPPORT_MPI)
  int numberOfNodes;
  MPI_Comm_size(MPI_COMM_WORLD, &numberOfNodes);

  if (numberOfNodes < 0) {
    throw std::runtime_error(
        "getNumberOfNode: Please check your MPI installation, as we got a "
        "negative number of nodes.");
  }

  return numberOfNodes;
#else
  return 1;
#endif
}

arma::uword nodeRank() {
#if defined(SUPPORT_MPI)
  int nodeRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeRank);

  if (nodeRank < 0) {
    throw std::runtime_error(
        "getNodeRank: Please check your MPI installation, as we got a negative "
        "node rank.");
  }

  return nodeRank;
#else
  return 0;
#endif
}

}  // namespace pass
