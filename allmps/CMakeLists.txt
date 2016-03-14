set(target MCAll)

set(srcs
  Algorithm/DuplicateRemover.cpp
  Implementations/MallAlgo.cpp
  )

find_package(OpenMP)

if(NOT OPENMP_FOUND)
  message(SEND_ERROR
    "Could not find a compatible OpenMP compiler. Consider turning BUILD_OPENMP to OFF")
endif()

find_package(MPI)

if(NOT MPI_CXX_FOUND)
  message(SEND_ERROR
    "Could not find a compatible MPI compiler. Consider turning BUILD_MPI to OFF")
endif()

add_executable(${target} main.cpp ${srcs})

target_link_libraries(${target}
  mantevo_marching_cubes_common
  ${MPI_CXX_LIBRARIES}
  )
target_compile_options(${target}
  PUBLIC ${OpenMP_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}
  )
set_target_options(${target}
  PROPERTIES LINK_FLAGS ${MPI_CXX_LINK_FLAGS}
  )
