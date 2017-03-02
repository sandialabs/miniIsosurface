cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 11)

project(flyingEdges)

find_package(CUDA REQUIRED)

add_subdirectory(serial)
add_subdirectory(cuda)
