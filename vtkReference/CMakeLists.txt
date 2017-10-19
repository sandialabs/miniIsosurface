# miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
# See LICENSE.txt for details.

# Copyright (c) 2017
# National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
# the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
# certain rights in this software.

cmake_minimum_required(VERSION 2.8)

PROJECT(vtkReference)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(srcs
    ../marchingCubes/util/Timer.cpp
    ../marchingCubes/mantevoCommon/YAML_Element.cpp
    ../marchingCubes/mantevoCommon/YAML_Doc.cpp
    )

add_executable(vtkReference MACOSX_BUNDLE main.cpp ${srcs})

if(VTK_LIBRARIES)
  target_link_libraries(vtkReference ${VTK_LIBRARIES})
else()
  target_link_libraries(vtkReference vtkHybrid vtkWidgets)
endif()
