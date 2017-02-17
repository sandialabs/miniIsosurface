set(target flyingEdgesSerial)

set(srcs
    FlyingEdgesAlgorithm.cpp
    ../util/Image3D.cpp
    ../util/Timer.cpp
    ../mantevoCommon/YAML_Doc.cpp
    ../mantevoCommon/YAML_Element.cpp
    )

add_executable(${target} main.cpp ${srcs})


