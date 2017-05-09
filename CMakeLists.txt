cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(PASS CXX)

find_path(MANTELLA_INCLUDE_DIR mantella0)

add_executable(example
  src/optimiser/rfc_1149_5.cpp
  src/main.cpp
)
set_property(TARGET example PROPERTY CXX_STANDARD 14)
set_property(TARGET example PROPERTY CXX_STANDARD_REQUIRED ON)
target_include_directories(example SYSTEM PUBLIC ${MANTELLA_INCLUDE_DIR})
