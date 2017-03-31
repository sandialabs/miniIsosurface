#ifndef FE_UTIL_STUFF
#define FE_UTIL_STUFF

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "Thrust_Config.h"

namespace util{

struct step
  : public thrust::unary_function<
        int,
        int>
{
    step(int const step_size)
      : step_size(step_size)
    {}

    __host__ __device__
    int
    operator()(int const& idx) const
    {
        return idx*step_size;
    }

    int const step_size;
};

template <typename Iter>
struct get_iter
  : public thrust::unary_function<
        int,
        Iter>
{
    get_iter(
        Iter iter)
      : iter(iter)
    {}

    __host__ __device__
    Iter
    operator()(int const& idx) const
    {
        return iter + idx;
    }

    Iter iter;
};

template <typename T>
using get_ptr = get_iter<typename thrust::device_vector<T>::const_pointer>;

}

#endif
