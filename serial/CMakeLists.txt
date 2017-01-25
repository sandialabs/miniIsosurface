set(target serial)

set(srcs
    ../util/Image3D.cpp
    )

add_executable(${target} main.cpp ${srcs})


