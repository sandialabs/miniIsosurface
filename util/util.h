#ifndef UTIL_UTIL_UTIL
#define UTIL_UTIL_UTIL

#include <array>

#include "MarchingCubesTables.h"

namespace util
{
    static const int caseMask[] = {1, 2, 4, 8, 16, 32, 64, 128};

    template<typename T>
    int
    findCaseId(std::array<T, 8> const& cubeVertexVals, T const& isoval)
    {
        int caseId=0;
        for(int i = 0; i < 8; ++i)
        {
            // Note that util::caseMask[i] = {1, 2, 4, 8, 16, 32, 64, 128}
            // Example:
            //   Suppose caseId = 3, i = 4 and cubeVertexVals[4] >= isoval.
            //   Then caseId will be set to caseId | 16. In terms of bits,
            //   caseId = 11000000 | 00001000 = 11001000 = 19
            if(cubeVertexVals[i] >= isoval)
            {
                caseId |= caseMask[i];
            }
        }
        return caseId;
    }

    template <typename T, std::size_t N>
    std::array<T, N>
    interpolate(std::array<T, N> const& a, std::array<T, N> const& b, T weight)
    {
        std::array<T, N> ret;
        for(int i = 0; i != N; ++i)
        {
            ret[i] = a[i] + (weight * (b[i] - a[i]));
        }
        return ret;
    }
}

#endif
