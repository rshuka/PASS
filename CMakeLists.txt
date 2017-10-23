# Require minimum Version of CMAKE
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

project(PASS LANGUAGES CXX)

find_package(Armadillo REQUIRED 7.800.2)

# -------------------------
# CONFIGURING LIBRARY BUILD
# -------------------------

message(STATUS "")
message(STATUS "Configuring library build.")

# All paths must start with "src/"

add_executable(pass
  src/helper/geometry.cpp
  src/helper/random.cpp
  src/helper/gtoc1/Astro_Functions.cpp
  src/helper/gtoc1/astro_helpers.cpp
  src/helper/gtoc1/constants.cpp
  src/helper/gtoc1/Lambert.cpp
  src/helper/gtoc1/PowSwingByInv.cpp
  src/helper/gtoc1/vector3d_helpers.cpp
  src/helper/gtoc1/ZeroFinder.cpp
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

  src/optimiser.cpp
  src/optimiser/hooke_jeeves_algorithm.cpp
  src/optimiser/particle_swarm_optimisation.cpp
  src/optimiser/random_search.cpp

  exe/main.cpp
)
set_property(TARGET pass PROPERTY CXX_STANDARD 14)
set_property(TARGET pass PROPERTY CXX_STANDARD_REQUIRED ON)
target_include_directories(pass SYSTEM PUBLIC ${ARMADILLO_INCLUDE_DIR})
target_link_libraries(pass PUBLIC ${ARMADILLO_LIBRARIES})
