cmake_minimum_required(VERSION 2.8.12)

project(mantevo-marching-cubes)

option(BUILD_OPENMP OFF)
option(BUILD_MPI OFF)

add_subdirectory(common)
add_subdirectory(dataGen)

add_subdirectory(serial)

if(BUILD_OPENMP)
  add_subdirectory(openmp)
  add_subdirectory(mergemp)
  if(BUILD_MPI)
    add_subdirectory(allmps)
  endif()
endif()

if(BUILD_MPI)
  add_subdirectory(mpi)
endif()
