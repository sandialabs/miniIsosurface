cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 11)

project(marchingCubes)

option(BUILD_OPENMP OFF)
option(BUILD_MPI OFF)

add_subdirectory(serial)

if(BUILD_OPENMP)
    add_subdirectory(openmp)
    add_subdirectory(openmpDupFree)
endif()

if(BUILD_MPI)
    add_subdirectory(mpi)
endif()

if(BUILD_OPENMP AND BUILD_MPI)
    add_subdirectory(openmpAndMpi)
endif()



