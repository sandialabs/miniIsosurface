/*
 * Image3D.cpp
 *
 *  Created on: Feb 17, 2017
 *      Author: dbourge
 */

#include "Image3D.h"

namespace util {

std::vector<scalar_t>::const_iterator
Image3D::getRowIter(size_t j, size_t k) const
{
    return data.cbegin() + nx*(k*ny + j);
}

scalarCube_t
Image3D::getValsCube(size_t i, size_t j, size_t k) const
{
    scalarCube_t vals;

    size_t idx = i + j*nx + k*nx*ny;

    vals[0] = getData(i, j, k);
    vals[1] = getData(i + 1, j, k);
    vals[2] = getData(i + 1, j + 1, k);
    vals[3] = getData(i, j + 1, k);
    vals[4] = getData(i, j, k + 1);
    vals[5] = getData(i + 1, j, k + 1);
    vals[6] = getData(i + 1, j + 1, k + 1);
    vals[7] = getData(i, j + 1, k + 1);

    return vals;
}

cube_t
Image3D::getPosCube(size_t i, size_t j, size_t k) const
{
    cube_t pos;

    scalar_t xpos = zeroPos[0] + i * spacing[0];
    scalar_t ypos = zeroPos[1] + j * spacing[1];
    scalar_t zpos = zeroPos[2] + k * spacing[2];

    pos[0][0] = xpos;
    pos[0][1] = ypos;
    pos[0][2] = zpos;

    pos[1][0] = xpos + spacing[0];
    pos[1][1] = ypos;
    pos[1][2] = zpos;

    pos[2][0] = xpos + spacing[0];
    pos[2][1] = ypos + spacing[1];
    pos[2][2] = zpos;

    pos[3][0] = xpos;
    pos[3][1] = ypos + spacing[1];
    pos[3][2] = zpos;

    pos[4][0] = xpos;
    pos[4][1] = ypos;
    pos[4][2] = zpos + spacing[2];

    pos[5][0] = xpos + spacing[0];
    pos[5][1] = ypos;
    pos[5][2] = zpos + spacing[2];

    pos[6][0] = xpos + spacing[0];
    pos[6][1] = ypos + spacing[1];
    pos[6][2] = zpos + spacing[2];

    pos[7][0] = xpos;
    pos[7][1] = ypos + spacing[1];
    pos[7][2] = zpos + spacing[2];

    return pos;
}

cube_t
Image3D::getGradCube(size_t i, size_t j, size_t k) const
{
    cube_t grad;

    grad[0] = computeGradient(i, j, k);
    grad[1] = computeGradient(i + 1, j, k);
    grad[2] = computeGradient(i + 1, j + 1, k);
    grad[3] = computeGradient(i, j + 1, k);
    grad[4] = computeGradient(i, j, k + 1);
    grad[5] = computeGradient(i + 1, j, k + 1);
    grad[6] = computeGradient(i + 1, j + 1, k + 1);
    grad[7] = computeGradient(i, j + 1, k + 1);

    return grad;
}

///////////////////////////////////////////////////////////////////////////////
// Private helper functions
//////////////////////////////////////////////////////////////////////////////

inline scalar_t
Image3D::getData(size_t i, size_t j, size_t k) const
{
    return data[k*nx*ny + j*nx + i];
}

std::array<scalar_t, 3>
Image3D::computeGradient(size_t i, size_t j, size_t k) const
{
    std::array<std::array<scalar_t, 2>, 3> x;
    std::array<scalar_t, 3> run;

    size_t dataIdx = i + j*nx + k*nx*ny;

    if (i == 0)
    {
        x[0][0] = data[dataIdx + 1];
        x[0][1] = data[dataIdx];
        run[0] = spacing[0];
    }
    else if (i == (nx - 1))
    {
        x[0][0] = data[dataIdx];
        x[0][1] = data[dataIdx - 1];
        run[0] = spacing[0];
    }
    else
    {
        x[0][0] = data[dataIdx + 1];
        x[0][1] = data[dataIdx - 1];
        run[0] = 2 * spacing[0];
    }

    if (j == 0)
    {
        x[1][0] = data[dataIdx + nx];
        x[1][1] = data[dataIdx];
        run[1] = spacing[1];
    }
    else if (j == (ny - 1))
    {
        x[1][0] = data[dataIdx];
        x[1][1] = data[dataIdx - nx];
        run[1] = spacing[1];
    }
    else
    {
        x[1][0] = data[dataIdx + nx];
        x[1][1] = data[dataIdx - ny];
        run[1] = 2 * spacing[1];
    }

    if (k == 0)
    {
        x[2][0] = data[dataIdx + nx*ny];
        x[2][1] = data[dataIdx];
        run[2] = spacing[2];
    }
    else if (k == (nz - 1))
    {
        x[2][0] = data[dataIdx];
        x[2][1] = data[dataIdx - nx*ny];
        run[2] = spacing[2];
    }
    else
    {
        x[2][0] = data[dataIdx + nx*ny];
        x[2][1] = data[dataIdx - nx*ny];
        run[2] = 2 * spacing[2];
    }

    std::array<scalar_t, 3> ret;

    ret[0] = (x[0][1] - x[0][0]) / run[0];
    ret[1] = (x[1][1] - x[1][0]) / run[1];
    ret[2] = (x[2][1] - x[2][0]) / run[2];

    return ret;
}

///////////////////////////////////////////////////////////////////////////////

}
