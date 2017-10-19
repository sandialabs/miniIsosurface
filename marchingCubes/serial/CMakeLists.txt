# miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
# See LICENSE.txt for details.

# Copyright (c) 2017
# National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
# the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
# certain rights in this software.

set(target serial)

set(srcs
    ../util/Image3D.cpp
    ../util/Timer.cpp
    ../mantevoCommon/YAML_Element.cpp
    ../mantevoCommon/YAML_Doc.cpp
    )

add_executable(${target} main.cpp ${srcs})


