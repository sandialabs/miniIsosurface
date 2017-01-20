/*
 * Image3D.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: dbourge, sjmunn
 */

#include "Image3D.h"

template<typename T>
std::array<std::array<T, 3>, 8>
Image3D<T>::getPosCube(unsigned xidx, unsigned yidx, unsigned zidx) const
{
    std::array<std::array<T, 3>, 8> pos;

    T xpos = origin[0] + xidx * spacing[0];
    T ypos = origin[1] + yidx * spacing[1];
    T zpos = origin[2] + zidx * spacing[2];

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

template <typename T>
std::array<std::array<T, 3>, 8>
Image3D<T>::getGradCube(unsigned xidx, unsigned yidx, unsigned zidx) const
{
    std::array<std::array<T, 3>, 8> grad;

    grad[0] = computeGradient(xidx, yidx, zidx);
    grad[1] = computeGradient(xidx + 1, yidx, zidx);
    grad[2] = computeGradient(xidx + 1, yidx + 1, zidx);
    grad[3] = computeGradient(xidx, yidx + 1, zidx);
    grad[4] = computeGradient(xidx, yidx, zidx + 1);
    grad[5] = computeGradient(xidx + 1, yidx, zidx + 1);
    grad[6] = computeGradient(xidx + 1, yidx + 1, zidx + 1);
    grad[7] = computeGradient(xidx, yidx + 1, zidx + 1);

    return grad;
}

template <typename T>
std::array<T, 3>
Image3D<T>::computeGradient(unsigned xidx, unsigned yidx, unsigned zidx) const
{
    std::array<std::array<T, 2>, 3> x;
    std::array<T, 3> run;

    unsigned dataIdx = xidx + yidx * dim[0] + zidx * dim[0] * dim[1];
    if (xidx == 0)
    {
        x[0][0] = data[dataIdx + 1];
        x[0][1] = data[dataIdx];
        run[0] = spacing[0];
    }
    else if (xidx == (dim[0] - 1))
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

    if (yidx == 0)
    {
        x[1][0] = data[dataIdx + dim[0]];
        x[1][1] = data[dataIdx];
        run[1] = spacing[1];
    }
    else if (yidx == (dim[1] - 1))
    {
        x[1][0] = data[dataIdx];
        x[1][1] = data[dataIdx - dim[0]];
        run[1] = spacing[1];
    }
    else
    {
        x[1][0] = data[dataIdx + dim[0]];
        x[1][1] = data[dataIdx - dim[0]];
        run[1] = 2 * spacing[1];
    }

    if (zidx == 0)
    {
        x[2][0] = data[dataIdx + dim[0]*dim[1]];
        x[2][1] = data[dataIdx];
        run[2] = spacing[2];
    }
    else if (zidx == (dim[2] - 1))
    {
        x[2][0] = data[dataIdx];
        x[2][1] = data[dataIdx - dim[0]*dim[1]];
        run[2] = spacing[2];
    }
    else
    {
        x[2][0] = data[dataIdx + dim[0]*dim[1]];
        x[2][1] = data[dataIdx - dim[0]*dim[1]];
        run[2] = 2 * spacing[2];
    }

    std::array<T, 3> ret;

    ret[0] = (x[0][1] - x[0][0]) / run[0];
    ret[1] = (x[1][1] - x[1][0]) / run[1];
    ret[2] = (x[2][1] - x[2][0]) / run[2];

    return ret;
}

template <typename T>
unsigned
Image3D<T>::getGlobalEdgeIndex(unsigned x, unsigned y, unsigned z,
                               unsigned cubeEdgeIdx) const
{
    unsigned edgeIndex = 0;
    switch (cubeEdgeIdx) {
    // Edges parallel to the x-axis
    case 0:
        edgeIndex = edgeIndexXaxis(x, y, z);
        break;
    case 2:
        edgeIndex = edgeIndexXaxis(x, y + 1, z);
        break;
    case 4:
        edgeIndex = edgeIndexXaxis(x, y, z + 1);
        break;
    case 6:
        edgeIndex = edgeIndexXaxis(x, y + 1, z + 1);
        break;
        // Edges parallel to the y-axis
    case 3:
        edgeIndex = edgeIndexYaxis(x, y, z);
        break;
    case 1:
        edgeIndex = edgeIndexYaxis(x + 1, y, z);
        break;
    case 7:
        edgeIndex = edgeIndexYaxis(x, y, z + 1);
        break;
    case 5:
        edgeIndex = edgeIndexYaxis(x + 1, y, z + 1);
        break;
        // Edges parallel to the z-axis
    case 8:
        edgeIndex = edgeIndexZaxis(x, y, z);
        break;
    case 9:
        edgeIndex = edgeIndexZaxis(x + 1, y, z);
        break;
    case 10:
        edgeIndex = edgeIndexZaxis(x, y + 1, z);
        break;
    case 11:
        edgeIndex = edgeIndexZaxis(x + 1, y + 1, z);
        break;
    }
    return edgeIndex;
}

template<typename T>
unsigned
Image3D<T>::edgeIndexXaxis(unsigned x, unsigned y, unsigned z) const
{
    //unsigned nXedges = (dim[0] - 1) * dim[1] * dim[2];
    //unsigned nYedges = dim[0] * (dim[1] - 1) * dim[2];
    //unsigned nZedges = dim[0] * dim[1] * (dim[2] - 1);
    //
    //unsigned nXYedges = nXedges + nYedges;
    //unsigned nAllEdges = nXedges + nYedges + nZedges;
    //
    //unsigned rangeX = dim[0] - 1;
    //unsigned rangeY = dim[1] - 1;
    //unsigned rangeZ = dim[2] - 1;

    unsigned index = x + ((dim[0] - 1) * y) + (dim[1] * (dim[0] - 1) * z);
    return index;
}

template<typename T>
unsigned
Image3D<T>::edgeIndexYaxis(unsigned x, unsigned y, unsigned z) const
{
    unsigned nXedges = (dim[0] - 1) * dim[1] * dim[2];
    unsigned index = nXedges + x + (dim[0] * y) + (dim[0] * (dim[1] - 1) * z);
    return index;
}

template<typename T>
unsigned
Image3D<T>::edgeIndexZaxis(unsigned x, unsigned y, unsigned z) const
{
    unsigned nXedges = (dim[0] - 1) * dim[1] * dim[2];
    unsigned nYedges = dim[0] * (dim[1] - 1) * dim[2];
    unsigned nXYedges = nXedges + nYedges;

    unsigned index = nXYedges + x + (dim[0] * y) + (dim[0] * dim[1] * z);
    return index;
}

template <typename T>
typename Image3D<T>::Image3DBuffer
Image3D<T>::createBuffer(unsigned xbeg, unsigned yidx, unsigned zidx) const
{
    using Iter = typename std::vector<T>::const_iterator;

    unsigned bufferIdx = xbeg + (yidx * dim[0]) + (zidx * dim[0] * dim[1]);
    Iter beg = data.begin();

    Iter x1buffer = beg + bufferIdx;
    Iter x2buffer = beg + bufferIdx + dim[0];
    Iter x3buffer = beg + bufferIdx + dim[0] * dim[1];
    Iter x4buffer = beg + bufferIdx + dim[0] + dim[0] * dim[1];

    return Image3D<T>::Image3DBuffer(x1buffer, x2buffer, x3buffer, x4buffer, xbeg);
}

template <typename T>
std::array<T, 8> Image3D<T>::Image3DBuffer::getCubeVertexValues(unsigned xidx) const
{
    std::array<T, 8> cubeVertexVals;

    cubeVertexVals[0] = x1buffer[xidx - xBeginIdx];
    cubeVertexVals[1] = x1buffer[xidx - xBeginIdx + 1];

    cubeVertexVals[2] = x2buffer[xidx - xBeginIdx + 1];
    cubeVertexVals[3] = x2buffer[xidx - xBeginIdx];

    cubeVertexVals[4] = x3buffer[xidx - xBeginIdx];
    cubeVertexVals[5] = x3buffer[xidx - xBeginIdx + 1];

    cubeVertexVals[6] = x4buffer[xidx - xBeginIdx + 1];
    cubeVertexVals[7] = x4buffer[xidx - xBeginIdx];

    return cubeVertexVals;
}

// Need this for explicit instantiation
template class Image3D<double>;
template class Image3D<float>;
