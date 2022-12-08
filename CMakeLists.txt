cmake_minimum_required(VERSION 3.24)
project(mpi_defferential_solver)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(mpi_defferential_solver
        main.cpp)
