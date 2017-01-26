/*
 * TriangleMesh.h
 *
 * Created on: Jan 13, 2017
 *     Author: dbourge
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
        typename std::vector<std::array<unsigned, 3> >::const_iterator;

    TriangleMesh()
    {}

    TriangleMesh(
        std::vector<std::array<T, 3> > points,
        std::vector<std::array<T, 3> > normals,
        std::vector<std::array<unsigned, 3> > indexTriangles)
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
    std::vector<std::array<unsigned, 3> > indexTriangles;
};

}

#endif
