set(target dataGen)

set(srcs
    open-simplex-noise.c
    )

add_executable(${target} main.cpp ${srcs})


