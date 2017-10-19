/*
 * TriangleMesh.h
 *
 * Created on: Jan 13, 2017
 *     Author: dbourge
 *
 * miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
 * See LICENSE.txt for details.
 *
 * Copyright (c) 2017
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
 * the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
 * certain rights in this software.
 */

#ifndef TRIANGLEMESH_H_
#define TRIANGLEMESH_H_

#include <vector>
#include <array>

namespace util {

template <typename T>
class TriangleMesh
{
public:
    using PointIterator =
        typename std::vector<std::array<T, 3> >::const_iterator;
    using NormalIterator =
        typename std::vector<std::array<T, 3> >::const_iterator;
    using TriangleIterator =
        typename std::vector<std::array<size_t, 3> >::const_iterator;

    TriangleMesh()
    {}

    TriangleMesh(
        std::vector<std::array<T, 3> > points,
        std::vector<std::array<T, 3> > normals,
        std::vector<std::array<size_t, 3> > indexTriangles)
      : points(points),
        normals(normals),
        indexTriangles(indexTriangles)
    {}

    std::size_t numberOfVertices() const
    {
        return points.size();
    }

    std::size_t numberOfTriangles() const
    {
        return indexTriangles.size();
    }

    PointIterator pointsBegin() const
    {
        return points.begin();
    }

    PointIterator pointsEnd() const
    {
        return points.end();
    }

    NormalIterator normalsBegin() const
    {
        return normals.begin();
    }

    NormalIterator normalsEnd() const
    {
        return normals.end();
    }

    TriangleIterator trianglesBegin() const
    {
        return indexTriangles.begin();
    }

    TriangleIterator trianglesEnd() const
    {
        return indexTriangles.end();
    }

private:
    std::vector<std::array<T, 3> > points;
    std::vector<std::array<T, 3> > normals;
    std::vector<std::array<size_t, 3> > indexTriangles;
};

}

#endif
