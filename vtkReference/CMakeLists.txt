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
