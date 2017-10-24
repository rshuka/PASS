# Require minimum Version of CMAKE
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# Retrieves the version numbers from "include/pass", avoiding possible inconsistency.
file(STRINGS "./include/pass" PASS_HEADER_CONTENTS REGEX "#define PASS_VERSION_[A-Z]+ ")
string(REGEX REPLACE ".*#define PASS_VERSION_MAJOR ([0-9]+).*" "\\1" PASS_VERSION_MAJOR "${PASS_HEADER_CONTENTS}")
string(REGEX REPLACE ".*#define PASS_VERSION_MINOR ([0-9]+).*" "\\1" PASS_VERSION_MINOR "${PASS_HEADER_CONTENTS}")
string(REGEX REPLACE ".*#define PASS_VERSION_PATCH ([0-9]+).*" "\\1" PASS_VERSION_PATCH "${PASS_HEADER_CONTENTS}")
set(PASS_VERSION "${PASS_VERSION_MAJOR}.${PASS_VERSION_MINOR}.${PASS_VERSION_PATCH}")

# Project name
project(PASS VERSION "${PASS_VERSION}" LANGUAGES CXX)


# Adding additional support
option(SUPPORT_SIMD "Add SSE2, SSE3, SSE4, AVX, etc. support" ON)
option(SUPPORT_OPENMP "Add OpenMP support" ON)
option(SUPPORT_MPI "Add MPI support" OFF)

if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
endif()

if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
endif()

if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
endif()



message(STATUS "")
message(STATUS "Building PASS (version: ${PASS_VERSION}).")

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "- Setting build to 'Release'.")
  message(STATUS "  - Use 'cmake ... -DCMAKE_BUILD_TYPE=[Release|Debug]' to specify this.")
  set(CMAKE_BUILD_TYPE "Release")
endif()

# --------------------
# FINDING DEPENDENCIES
# --------------------

message(STATUS "")
message(STATUS "Finding dependencies.")

# Requite minimum Version of Amardillo
set(REQUIRED_ARMADILLO_VERSION "8.200.0")
find_package(Armadillo QUIET)
include_directories(${ARMADILLO_INCLUDE_DIR})
message(STATUS "Found Armadillo: ${ARMADILLO_LIBRARIES} (found suitable version ${ARMADILLO_VERSION_STRING}, minimum required is ${REQUIRED_ARMADILLO_VERSION}).")

# Support MPI
if (SUPPORT_MPI)
  message(STATUS "")
  message(STATUS "- Adding MPI support.")
  message(STATUS "  - Use 'cmake ... -DSUPPORT_MPI=Off' to exclude this.")
  find_package(MPI REQUIRED)
else()
  message(STATUS "- Excluding MPI support.")
  message(STATUS "  - Use 'cmake ... -DSUPPORT_MPI=ON' to add this.")
endif()

# Support OPENMP
if (SUPPORT_OPENMP)
  message(STATUS "")
  message(STATUS "- Trying to add OpenMP support.")
  message(STATUS "  - Use 'cmake ... -DSUPPORT_OPENMP=Off' to exclude this.")

  # Makes sure that `MAXIMAL_NUMBER_OF_THREADS` contains a sane value.
  if (NOT MAXIMAL_NUMBER_OF_THREADS OR NOT MAXIMAL_NUMBER_OF_THREADS MATCHES "^[1-9][0-9]*$")
    message(STATUS "  - Trying to determine the number of cores.")

    include(ProcessorCount)
    ProcessorCount(MAXIMAL_NUMBER_OF_THREADS)
    message(STATUS "    - Found ${MAXIMAL_NUMBER_OF_THREADS} core(s).")
  endif()

  if (MAXIMAL_NUMBER_OF_THREADS EQUAL 0)
    message(STATUS "- Excluding OpenMP support as the maximal number of threads to support was either not automatically determined or manually set to an integer less than 1.")
    message(STATUS "  - Use 'cmake ... -DMAXIMAL_NUMBER_OF_THREADS=...' to specify the maximal number of threads (>= 1) to support manually.")
    set(SUPPORT_OPENMP OFF)
  else()
    find_package(OpenMP)
    if (OPENMP_FOUND)
      message(STATUS "  - Adding OpenMP support (optimising for ${MAXIMAL_NUMBER_OF_THREADS} threads(s)).")
    message(STATUS "  - Use 'cmake ... -DMAXIMAL_NUMBER_OF_THREADS=...' to specify the maximal number of threads (>= 1) manually.")
      message(STATUS "  - Alternatively, set `#define MAXIMAL_NUMBER_OF_THREADS ...` in your source code, prior to including PASS the first time.")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else()
      message(STATUS "  - OpenMP is not supported by your compiler.")
      set(SUPPORT_OPENMP OFF)
    endif()
  endif()
else()
  message(STATUS "- Excluding OpenMP support.")
  message(STATUS "  - Use 'cmake ... -DSUPPORT_OPENMP=ON' to add this.")
endif()

# Makes sure that `MAXIMAL_NUMBER_OF_THREADS` contains a sane value, even if OpenMP is not to be supported.
if (NOT SUPPORT_OPENMP)
  set(MAXIMAL_NUMBER_OF_THREADS 1)
endif()

# -------------------------
# CONFIGURING LIBRARY BUILD
# -------------------------

message(STATUS "")
message(STATUS "Configuring library build.")

configure_file(${PROJECT_SOURCE_DIR}/include/pass_bits/config.hpp.cmake ${PROJECT_SOURCE_DIR}/include/pass_bits/config.hpp)

# All paths must start with "src/"

add_library(pass SHARED
  # Configuration
  src/config.cpp

  # Helper
  src/helper/geometry.cpp
  src/helper/random.cpp
  src/helper/gtoc1/Astro_Functions.cpp
  src/helper/gtoc1/astro_helpers.cpp
  src/helper/gtoc1/constants.cpp
  src/helper/gtoc1/Lambert.cpp
  src/helper/gtoc1/PowSwingByInv.cpp
  src/helper/gtoc1/vector3d_helpers.cpp
  src/helper/gtoc1/ZeroFinder.cpp

  # Optimisation Problems
  src/problem.cpp
  src/problem/ackley_function.cpp
  src/problem/gtoc1.cpp
  src/problem/parallel_kinematic_machine_3PRPR.cpp
  src/problem/parallel_kinematic_machine_3PRRR.cpp
  src/problem/parallel_kinematic_machine_6PRUS.cpp
  src/problem/parallel_kinematic_machine_6PUPS.cpp
  src/problem/rastrigin_function.cpp
  src/problem/rosenbrock_function.cpp
  src/problem/sphere_function.cpp
  src/problem/sum_of_different_powers.cpp

  # Optimisation algorithms
  src/optimiser.cpp
  src/optimiser/hooke_jeeves_algorithm.cpp
  src/optimiser/particle_swarm_optimisation.cpp
  src/optimiser/random_search.cpp
  src/optimiser/parallel_swarm_search/parallel_swarm_search.cpp
)
set_property(TARGET pass PROPERTY CXX_STANDARD 14)
set_property(TARGET pass PROPERTY CXX_STANDARD_REQUIRED ON)

if (SUPPORT_SIMD)
  message(STATUS "- Adding SSE3, SSE4, AVX, ... support.")
  message(STATUS "  - Use 'cmake ... -DSUPPORT_SIMD=Off' to exclude this.")
  target_compile_options(pass PRIVATE -march=native)
else()
  message(STATUS "- Excluding SSE3, SSE4, AVX, ... support.")
  message(STATUS "  - Use 'cmake ... -DSUPPORT_SIMD=ON' to add this.")
endif()

message(STATUS "- Using level 3 code optimisation.")
target_compile_options(pass PRIVATE -O3)

# Include PASS

target_include_directories(pass PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

target_include_directories(pass SYSTEM PUBLIC ${ARMADILLO_INCLUDE_DIR})
target_link_libraries(pass PUBLIC ${ARMADILLO_LIBRARIES})

if (SUPPORT_MPI)
  target_include_directories(pass SYSTEM PUBLIC ${MPI_INCLUDE_PATH})
  target_link_libraries(pass PUBLIC ${MPI_C_LIBRARIES})
endif()

# -------------
# BUILD LIBRARY
# -------------

message(STATUS "")
message(STATUS "Configuring library installation.")

if (NOT INSTALL_INCLUDE_DIR)
  set(INSTALL_INCLUDE_DIR "include")
endif()

set_property(TARGET pass PROPERTY VERSION ${PASS_VERSION})
set_property(TARGET pass PROPERTY SOVERSION ${PASS_VERSION_MAJOR})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev)

if (NOT INSTALL_LIB_DIR)
  set(INSTALL_LIB_DIR "lib")
endif()

if (NOT INSTALL_BIN_DIR)
  set(INSTALL_BIN_DIR "bin")
endif()

install(TARGETS pass
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${INSTALL_BIN_DIR})

# ----------------
# Status Updates
# ----------------

message(STATUS "")
message(STATUS "Noticable variables:")
message(STATUS "- CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "- CMAKE_LIBRARY_OUTPUT_DIRECTORY = ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(STATUS "- CMAKE_ARCHIVE_OUTPUT_DIRECTORY = ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
message(STATUS "- CMAKE_RUNTIME_OUTPUT_DIRECTORY = ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "- CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
message(STATUS "- SUPPORT_SIMD = ${SUPPORT_SIMD}")
message(STATUS "- SUPPORT_OPENMP = ${SUPPORT_OPENMP}")
message(STATUS "- SUPPORT_MPI = ${SUPPORT_MPI}")
message(STATUS "- ARMADILLO_LIBRARIES = ${ARMADILLO_LIBRARIES}")
message(STATUS "- INSTALL_INCLUDE_DIR = ${INSTALL_INCLUDE_DIR}")
message(STATUS "- INSTALL_LIB_DIR = ${INSTALL_LIB_DIR}")
message(STATUS "- INSTALL_BIN_DIR = ${INSTALL_BIN_DIR}")
message(STATUS "")
