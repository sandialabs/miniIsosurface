set(target serial)

set(srcs
    ../util/Image3D.cpp
    ../util/Timer.cpp
    ../mantevoCommon/YAML_Element.cpp
    ../mantevoCommon/YAML_Doc.cpp
    )

add_executable(${target} main.cpp ${srcs})


