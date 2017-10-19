# miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
# See LICENSE.txt for details.

# Copyright (c) 2017
# National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
# the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
# certain rights in this software.

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



