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

using std::size_t;

namespace util {

template <typename T>
class Image3D
{
public:

    // This constructor is used to construct images that contain
    // only a section of data. When this is the case, data must contain
    // more points than just indexed by indexBeg to indexEnd. This is because
    // calling getPosCube or getGradCube require values outside of
    // indexBeg and indexEnd.
    Image3D(std::vector<T> data,
            std::array<T, 3> spacing,
            std::array<T, 3> zeroPos,
            std::array<size_t, 3> indexBeg,
            std::array<size_t, 3> indexEnd,
            std::array<size_t, 3> dataBeg,
            std::array<size_t, 3> dataEnd,
            std::array<size_t, 3> globalDim)
      : data(data), spacing(spacing), zeroPos(zeroPos),
        indexBeg(indexBeg), indexEnd(indexEnd),
        dataBeg(dataBeg), dataEnd(dataEnd),
        globalDim(globalDim)
    {}

    // This constructor is used to construct an image of size
    // dimensions.
    Image3D(std::vector<T> data,
            std::array<T, 3> spacing,
            std::array<T, 3> zeroPos,
            std::array<size_t, 3> dimensions)
      : data(data), spacing(spacing), zeroPos(zeroPos),
        indexBeg({0, 0, 0}),
        indexEnd({dimensions[0] - 1, dimensions[1] - 1, dimensions[2] - 1}),
        dataBeg({0, 0, 0}), dataEnd(dimensions),
        globalDim(dimensions)
    {}

    std::array<std::array<T, 3>, 8>
    getPosCube(size_t xidx, size_t yidx, size_t zidx) const;

    std::array<std::array<T, 3>, 8>
    getGradCube(size_t xidx, size_t yidx, size_t zidx) const;

    size_t
    getGlobalEdgeIndex(size_t xidx, size_t yidx, size_t zidx,
                       size_t cubeEdgeIdx) const;

    size_t xBeginIdx() const { return indexBeg[0]; }
    size_t yBeginIdx() const { return indexBeg[1]; }
    size_t zBeginIdx() const { return indexBeg[2]; }
    size_t xEndIdx()   const { return indexEnd[0]; }
    size_t yEndIdx()   const { return indexEnd[1]; }
    size_t zEndIdx()   const { return indexEnd[2]; }

    size_t xdimension() const { return globalDim[0]; }
    size_t ydimension() const { return globalDim[1]; }
    size_t zdimension() const { return globalDim[2]; }

    // TODO
    void cut_down(int const& i)
    {
        data.resize(i*ydimension()*zdimension());
        dataEnd[0] = i;
        globalDim[0] = i;
        indexEnd[0] = i;
    }

private:
    std::array<T, 3>
    computeGradient(size_t xidx, size_t yidx, size_t zidx) const;

    size_t
    edgeIndexXaxis(size_t x, size_t y, size_t z) const;

    size_t
    edgeIndexYaxis(size_t x, size_t y, size_t z) const;

    size_t
    edgeIndexZaxis(size_t x, size_t y, size_t z) const;

    class Image3DBuffer
    {
    private:
        using Iter = typename std::vector<T>::const_iterator;
    public:
        Image3DBuffer(
            Iter x1buffer, Iter x2buffer,
            Iter x3buffer, Iter x4buffer,
            size_t xBeginIdx)
          : x1buffer(x1buffer), x2buffer(x2buffer),
            x3buffer(x3buffer), x4buffer(x4buffer),
            xBeginIdx(xBeginIdx)
        {}

        std::array<T, 8> getCubeVertexValues(size_t xidx) const;

    private:
        Iter x1buffer;
        Iter x2buffer;
        Iter x3buffer;
        Iter x4buffer;
        size_t const xBeginIdx;
    };
public:
    Image3DBuffer
    createBuffer(size_t xbeg, size_t yidx, size_t zidx) const;

private:
    std::vector<T>          data;       // A vector containing scalar values
                                        // along three-dimensional space.

    std::array<T, 3>        spacing;    // The distance between two points in
                                        // the mesh.

    std::array<T, 3>        zeroPos;    // The position at index (0, 0, 0).

    std::array<size_t, 3> indexBeg;     // The indices that can be used to get
    std::array<size_t, 3> indexEnd;     // information for a cube are all indices
                                        // such that index[i] is in
                                        // [indexBeg[i], indexEnd[i]) where i
                                        // is 0, 1, or 2.

    std::array<size_t, 3> dataBeg;      // The indices that data contains are
    std::array<size_t, 3> dataEnd;      // given by dataBeg and dataEnd. So
                                        // all indices such that index[i] is in
                                        // [dataBeg[i], dataEnd[i]) where i
                                        // is 0, 1 or 2.

    std::array<size_t, 3> globalDim;    // The dimension of the entire
                                        // image.
};

} // util namespace

#endif
