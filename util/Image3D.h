/*
 * Image3D.h
 *
 *  Created on: Jan 13, 2017
 *      Author: dbourge, sjmunn
 */

#ifndef IMAGE3DREADER_H_
#define IMAGE3DREADER_H_

#include <array>
#include <vector>

//#include <iostream>

namespace util {

template <typename T>
class Image3D
{
public:
    Image3D(std::vector<T> data,
            std::array<unsigned, 3> dim,
            std::array<T, 3> spacing,
            std::array<T, 3> origin)
        : data(data), dim(dim),
        spacing(spacing), origin(origin),
        indexOrigin({0, 0, 0})
    {}

    Image3D(std::vector<T> data,
            std::array<unsigned, 3> dim,
            std::array<T, 3> spacing,
            std::array<T, 3> origin,
            std::array<unsigned, 3> indexOrigin)
        : data(data), dim(dim),
        spacing(spacing), origin(origin),
        indexOrigin(indexOrigin)
    {}

    std::array<std::array<T, 3>, 8>
    getPosCube(unsigned xidx, unsigned yidx, unsigned zidx) const;

    std::array<std::array<T, 3>, 8>
    getGradCube(unsigned xidx, unsigned yidx, unsigned zidx) const;

    unsigned
    getGlobalEdgeIndex(unsigned xidx, unsigned yidx, unsigned zidx,
                 unsigned cubeEdgeIdx) const;

    unsigned xBeginIdx() const { return indexOrigin[0]; }
    unsigned yBeginIdx() const { return indexOrigin[1]; }
    unsigned zBeginIdx() const { return indexOrigin[2]; }
    unsigned xEndIdx() const { return indexOrigin[0] + dim[0]; }
    unsigned yEndIdx() const { return indexOrigin[1] + dim[1]; }
    unsigned zEndIdx() const { return indexOrigin[2] + dim[2]; }

private:
    std::array<T, 3>
    computeGradient(unsigned xidx, unsigned yidx, unsigned zidx) const;

    unsigned
    edgeIndexXaxis(unsigned x, unsigned y, unsigned z) const;

    unsigned
    edgeIndexYaxis(unsigned x, unsigned y, unsigned z) const;

    unsigned
    edgeIndexZaxis(unsigned x, unsigned y, unsigned z) const;

    class Image3DBuffer
    {
    private:
        using Iter = typename std::vector<T>::const_iterator;
    public:
        Image3DBuffer(
            Iter x1buffer, Iter x2buffer,
            Iter x3buffer, Iter x4buffer,
            unsigned xBeginIdx)
          : x1buffer(x1buffer), x2buffer(x2buffer),
            x3buffer(x3buffer), x4buffer(x4buffer),
            xBeginIdx(xBeginIdx)
        {}

        std::array<T, 8> getCubeVertexValues(unsigned xidx) const;

    private:
        Iter x1buffer;
        Iter x2buffer;
        Iter x3buffer;
        Iter x4buffer;
        unsigned const xBeginIdx;
    };
public:
    Image3DBuffer
    createBuffer(unsigned xbeg, unsigned yidx, unsigned zidx) const;

private:
    std::vector<T>          data;
    std::array<unsigned, 3> dim;
    std::array<T, 3>        spacing;
    std::array<T, 3>        origin;
    std::array<unsigned, 3> indexOrigin;
};

} // util namespace

#endif
