cmake_minimum_required(VERSION 3.0.0)
project(fluidic_simulations)

include(CTest)
enable_testing()

add_executable(fluidic_simulations main.cpp Simulation/simulation.c)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

