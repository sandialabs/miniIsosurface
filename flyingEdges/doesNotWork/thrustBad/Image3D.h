#ifndef IMAGE3D_H
#define IMAGE3D_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "Thrust_Config.h"

struct Image3D
{
public:
    Image3D(thrust::host_vector<scalar_t> data_,
            scalar_t const& spacingX,
            scalar_t const& spacingY,
            scalar_t const& spacingZ,
            scalar_t const& zeroPosX,
            scalar_t const& zeroPosY,
            scalar_t const& zeroPosZ,
            int const& nx,
            int const& ny,
            int const& nz)
      : data(data_),
        spacing(spacingX, spacingY, spacingZ),
        zeroPos(zeroPosX, zeroPosY, zeroPosZ),
        n(nx, ny, nz)
    {}

    int xdimension() const { return n.x; }
    int ydimension() const { return n.y; }
    int zdimension() const { return n.z; }

    using iterator =
        typename thrust::device_vector<scalar_t>::const_iterator;

    iterator begin() const
    {
        return data.begin();
    }

    iterator end() const
    {
        return data.end();
    }

    thrust::device_vector<scalar_t> const data;
    t3<const scalar_t> spacing;
    t3<const scalar_t> zeroPos;
    t3<const int> n;
};

#endif
