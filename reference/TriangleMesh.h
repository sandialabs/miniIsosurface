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

template <typename T>
class TriangleMesh
{
public:
    TriangleMesh(
        std::vector<std::array<T, 3> > points,
        std::vector<std::array<T, 3> > normals,
        std::vector<std::array<unsigned, 3> > indexTriangles)
      : points(points),
        normals(normals),
        indexTriangles(indexTriangles)
    {}

private:
    std::vector<std::array<T, 3> > points;
    std::vector<std::array<T, 3> > normals;
    std::vector<std::array<unsigned, 3> > indexTriangles;
};

#endif
