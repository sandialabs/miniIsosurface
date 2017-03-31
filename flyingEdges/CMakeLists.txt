cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 11)

project(flyingEdges)

find_package(CUDA REQUIRED)

add_subdirectory(serial)
add_subdirectory(thrust)
#add_subdirectory(cuda)
#add_subdirectory(cuda2)
#add_subdirectory(exp)

