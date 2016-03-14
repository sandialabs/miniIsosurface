set(target MCSerial)

set(srcs
  Implementations/SerialAlgo.cpp
  )

add_executable(${target} main.cpp ${srcs})

target_link_libraries(${target} manteveo_marching_cubes_common)
