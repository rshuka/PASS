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
  src/helper/random.cpp
  src/problem.cpp
  src/problem/ackley_function.cpp
  src/problem/sphere_function.cpp

  src/optimiser.cpp
  src/optimiser/particle_swarm_optimisation.cpp
  src/optimiser/random_search.cpp

  exe/main.cpp
)
set_property(TARGET pass PROPERTY CXX_STANDARD 14)
set_property(TARGET pass PROPERTY CXX_STANDARD_REQUIRED ON)
target_include_directories(pass SYSTEM PUBLIC ${ARMADILLO_INCLUDE_DIR})
target_link_libraries(pass PUBLIC ${ARMADILLO_LIBRARIES})
