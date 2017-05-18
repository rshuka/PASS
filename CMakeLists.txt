cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(PASS CXX)

find_path(MANTELLA_INCLUDE_DIR mantella0)
find_package(Armadillo REQUIRED 7.300.0)

add_executable(pass
  src/problem.cpp
  src/problem/sphere_function.cpp

  src/optimiser.cpp
  src/optimiser/random_search.cpp

  src/main.cpp
)
set_property(TARGET pass PROPERTY CXX_STANDARD 14)
set_property(TARGET pass PROPERTY CXX_STANDARD_REQUIRED ON)
target_include_directories(pass SYSTEM PUBLIC ${ARMADILLO_INCLUDE_DIR})
target_include_directories(pass SYSTEM PUBLIC ${MANTELLA_INCLUDE_DIR})
target_link_libraries(pass PUBLIC ${ARMADILLO_LIBRARIES})
