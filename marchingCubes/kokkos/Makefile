# miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
# See LICENSE.txt for details.

# Copyright (c) 2017
# National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
# the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
# certain rights in this software.

KOKKOS_PATH = ${HOME}/kokkos
KOKKOS_DEVICES = "OpenMP"
EXE_NAME = "kokkos"

SRC = $(wildcard *.cpp)

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/config/nvcc_wrapper
EXE = ${EXE_NAME}.cuda
KOKKOS_ARCH = "SNB,Kepler35"
KOKKOS_CUDA_OPTIONS = "enable_lambda"
else
CXX = g++
EXE = ${EXE_NAME}.host
KOKKOS_ARCH = "SNB"
endif

CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS =

DEPFLAGS = -M

OBJ = $(SRC:.cpp=.o) Image3D.o Timer.o YAML_Doc.o YAML_Element.o
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(MMC_OBJ) $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(MMC_OBJ) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

Image3D.o : ../util/Image3D.h ../util/Image3D.cpp
	$(CXX) $(CXXFLAGS) -c ../util/Image3D.cpp
#Image3D.o : Image3D.h Image3D.cpp
#	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c Image3D.cpp
Timer.o : ../util/Timer.h ../util/Timer.cpp
	$(CXX) $(CXXFLAGS) -c ../util/Timer.cpp
YAML_Doc.o : ../mantevoCommon/YAML_Doc.hpp ../mantevoCommon/YAML_Doc.cpp
	$(CXX) $(CXXFLAGS) -c ../mantevoCommon/YAML_Doc.cpp
YAML_Element.o : ../mantevoCommon/YAML_Element.hpp ../mantevoCommon/YAML_Element.cpp
	$(CXX) $(CXXFLAGS) -c ../mantevoCommon/YAML_Element.cpp
