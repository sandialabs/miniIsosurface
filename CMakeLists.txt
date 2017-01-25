cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 11)

project(mantevo-marching-cubes)

add_subdirectory(tests)
add_subdirectory(reference)
add_subdirectory(openmp)
add_subdirectory(openmpDupFree)
add_subdirectory(mpi)



