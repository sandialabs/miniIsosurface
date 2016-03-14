set(target MCoMP)

set(srcs
  Implementations/OpenMPAlgo.cpp
  )

find_package(OpenMP)

if (NOT OPENMP_FOUND)
  message(SEND_ERROR
    "Could not find a compatible OpenMP compiler. Consider turning BUILD_OPENMP to OFF")
endif()

add_executable(${target} main.cpp ${srcs})

target_link_libraries(${target} mantevo_marching_cubes_common)
target_compile_options(${target} PUBLIC ${OpenMP_CXX_FLAGS})
set_target_properties(${target} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
