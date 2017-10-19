# miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
# See LICENSE.txt for details.

# Copyright (c) 2017
# National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
# the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
# certain rights in this software.

set(target openmpAndMpi)

set(srcs
    ../util/Image3D.cpp
    ../util/Timer.cpp
    ../mantevoCommon/YAML_Element.cpp
    ../mantevoCommon/YAML_Doc.cpp
    )

find_package(MPI)

if (NOT MPI_CXX_FOUND)
  message(SEND_ERROR
      "Could not find a compatible MPI compiler. Consider turning BUILD_MPI to OFF")
endif()

find_package(OpenMP)

if (NOT OPENMP_FOUND)
  message(SEND_ERROR
      "Could not find a compatible OpenMP compiler. Consider turning BUILD_OPENMP to OFF")
endif()


add_executable(${target} main.cpp ${srcs})

target_include_directories(${target}
  PUBLIC "${MPI_CXX_INCLUDE_PATH}"
  )
target_compile_options(${target}
  PUBLIC ${OpenMP_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}
  )
target_link_libraries(${target}
  ${MPI_CXX_LIBRARIES}
  )
set_target_properties(${target}
  PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS} ${MPI_CXX_LINK_FLAGS}"
  )
