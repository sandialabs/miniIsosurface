cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 11)

project(flyingEdges)

option(BUILD_CUDA OFF)

add_subdirectory(serial)

if (BUILD_CUDA)
    find_package(CUDA REQUIRED)
    add_subdirectory(thrust)
    #add_subdirectory(cuda)
    #add_subdirectory(cuda2)
endif()
