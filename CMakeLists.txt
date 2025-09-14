cmake_minimum_required(VERSION 3.12.0)
project(project LANGUAGES CXX C)

include(${CMAKE_CURRENT_LIST_DIR}/../common/cmake-common.cmake)

add_executable(project src/main.cpp src/Poisson3D.cpp)
deal_ii_setup_target(project)