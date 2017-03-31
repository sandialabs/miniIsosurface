#ifndef CONFIG_STUFF
#define CONFIG_STUFF

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/scan.h>

#include <thrust/functional.h>

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>

#include <thrust/execution_policy.h>

#define RUN_ON_DEVICE

///////////////////////////////////////////////////////////////////////////////
#ifdef RUN_ON_DEVICE
#define __platform__ __device__
#else
#define __platform__ __host__
#endif

///////////////////////////////////////////////////////////////////////////////
// Basic typedef stuff
using thrust::host_vector;
using thrust::device_vector;

template <typename T>
#ifdef RUN_ON_DEVICE
using vector = thrust::device_vector<T>;
#else
using vector = thrust::host_vector<T>;
#endif

template <typename T>
using pointer = typename vector<T>::pointer;

template <typename T>
using const_pointer = typename vector<T>::const_pointer;

using scalar_t = double;
using uchar = unsigned char;

#ifdef RUN_ON_DEVICE
auto& policy = thrust::device;
#else
auto& policy = thrust::host;
#endif

///////////////////////////////////////////////////////////////////////////////
// Table stuff
#ifdef RUN_ON_DEVICE
#include "CudaMarchingCubesTables.h"
#define mctable cuda_util
#else
#include "MarchingCubesTables.h"
#define mctable util
#endif

///////////////////////////////////////////////////////////////////////////////
// Iterator stuff
using thrust::make_counting_iterator;
using thrust::make_zip_iterator;
using thrust::make_tuple;

///////////////////////////////////////////////////////////////////////////////
// Algorithm stuff
using thrust::transform;
using thrust::for_each;

///////////////////////////////////////////////////////////////////////////////
// Other thrust stuff
using thrust::unary_function;
using thrust::binary_function;
using thrust::tuple;
using thrust::get;

#endif
