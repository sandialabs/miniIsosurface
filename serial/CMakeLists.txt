set(target serial)

set(srcs
    ../util/Image3D.cpp
    ../util/Timer.cpp
    ../mantevo/YAML_Element.cpp
    ../mantevo/YAML_Doc.cpp
    )

add_executable(${target} main.cpp ${srcs})


