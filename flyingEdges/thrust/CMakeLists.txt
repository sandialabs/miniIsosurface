set(target flyingEdgesThrust)

set(srcs
    ../util/Timer.cpp
    ../mantevoCommon/YAML_Doc.cpp
    ../mantevoCommon/YAML_Element.cpp
    )

set(CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS} -std=c++11
    -O3
    )

cuda_add_executable(${target} main.cu ${srcs})
#add_executable(${target} main.cpp ${srcs})


