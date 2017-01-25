set(target mpi)

set(srcs
    ../util/Image3D.cpp
    )

find_package(MPI)

if (NOT MPI_CXX_FOUND)
  message(SEND_ERROR
    "Could not find a compatible MPI compiler.")
endif()

add_executable(${target} main.cpp ${srcs})

target_include_directories(${target}
  PUBLIC "${MPI_CXX_INCLUDE_PATH}"
  )
target_compile_options(${target}
  PUBLIC ${MPI_CXX_COMPILE_FLAGS}
  )
target_link_libraries(${target}
    ${MPI_CXX_LIBRARIES}
    )
set_target_properties(${target}
  PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}"
  )
