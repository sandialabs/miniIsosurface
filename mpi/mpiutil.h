/*
 * mpi/mpiutil.h
 *
 *  Created on: Feb 2, 2017
 *      Author: dbourge
 */

#ifndef MPIUTHELPERFUNCTIONSM
#define MPIUTHELPERFUNCTIONSM

#include <array>
#include <vector>
#include <unordered_map>

#include <algorithm>

#include "../util/TriangleMesh.h"

#include <mpi.h>

namespace mpiutil {

template<typename T, std::size_t N>
struct arrayHash
{
    std::size_t operator()(std::array<T, N> const &arr) const
    {
        std::size_t sum(0);
        for(auto &&i : arr)
            sum += std::hash<T>()(i);
        return sum;
    }
};

template <typename T>
using UnorderedMapArr =
    std::unordered_map<std::array<T, 3>, size_t, arrayHash<T, 3> >;

template <typename T>
util::TriangleMesh<T>
mergeMeshes(std::vector<util::TriangleMesh<T> > const& meshes)
{
    if(meshes.size() == 0)
    {
        return util::TriangleMesh<T>();
    }

    std::vector<std::array<T, 3> > points;
    std::vector<std::array<T, 3> > normals;
    std::vector<std::array<size_t, 3> > indexTriangles;

    UnorderedMapArr<T> pointMap;

    size_t count = 0;
    for(util::TriangleMesh<T> const& mesh: meshes)
    {
        // add to points and normals, filling up pointMap
        auto thesePoints = mesh.pointsBegin();
        auto theseNormals = mesh.normalsBegin();
        size_t numPoints = mesh.numberOfVertices();
        for(size_t idx = 0; idx != numPoints; ++idx)
        {
            std::array<T, 3> const& currentPoint = thesePoints[idx];
            std::array<T, 3> const& currentNormal = theseNormals[idx];

            if(pointMap.find(currentPoint) == pointMap.end())
            {
                pointMap[currentPoint] = count++;

                points.push_back(currentPoint);
                normals.push_back(currentNormal);
            }
        }

        // Use pointMap to transform the current mesh's indexTriangles
        // for the current global indexTriangles.
        auto triBeg = mesh.trianglesBegin();
        auto triEnd = mesh.trianglesEnd();
        for(auto triIt = triBeg; triIt != triEnd; ++triIt)
        {
            std::array<size_t, 3> const& oldTri = *triIt;

            std::array<size_t, 3> tri;
            tri[0] = pointMap[thesePoints[oldTri[0]]];
            tri[1] = pointMap[thesePoints[oldTri[1]]];
            tri[2] = pointMap[thesePoints[oldTri[2]]];

            indexTriangles.push_back(tri);
        }
    }

    return util::TriangleMesh<T>(points, normals, indexTriangles);
}

std::vector<util::TriangleMesh<float> >
gatherMeshes(util::TriangleMesh<float> const& thisMesh,
    int pid,
    std::vector<size_t> numVerts,
    std::vector<size_t> numTris)
{
    // Flatten this points, normals
    size_t thisNumVerts = thisMesh.numberOfVertices();
    auto thisPoints = thisMesh.pointsBegin();
    auto thisNormals = thisMesh.normalsBegin();

    std::vector<float> flatThisMesh(thisNumVerts * 6);
    for(size_t idx = 0; idx != thisNumVerts; ++idx)
    {
        flatThisMesh[idx * 6 + 0] = thisPoints[idx][0];
        flatThisMesh[idx * 6 + 1] = thisPoints[idx][1];
        flatThisMesh[idx * 6 + 2] = thisPoints[idx][2];
        flatThisMesh[idx * 6 + 3] = thisNormals[idx][0];
        flatThisMesh[idx * 6 + 4] = thisNormals[idx][1];
        flatThisMesh[idx * 6 + 5] = thisNormals[idx][2];
    }

    // Flatten this indexTriangles
    size_t thisNumTris = thisMesh.numberOfTriangles();
    auto thisTris = thisMesh.trianglesBegin();

    std::vector<size_t> flatThisTris(thisNumTris * 3);
    for(size_t idx = 0; idx != thisNumTris; ++idx)
    {
        flatThisTris[idx * 3 + 0] = thisTris[idx][0];
        flatThisTris[idx * 3 + 1] = thisTris[idx][1];
        flatThisTris[idx * 3 + 2] = thisTris[idx][2];
    }

    // Communicate data

    size_t numProcesses = numVerts.size();

    numVerts.push_back(0);
    numTris.push_back(0);
    size_t maxVerts = *std::max_element(numVerts.begin(), numVerts.end());
    size_t maxTris = *std::max_element(numTris.begin(), numTris.end());

    size_t receiveVertsSize = 3*2*maxVerts;
    size_t receiveTrisSize = 3*maxTris;

    std::vector<float> receiveVerts(numProcesses * receiveVertsSize);
    std::vector<size_t> receiveTris(numProcesses * receiveTrisSize);

    MPI_Gather(flatThisMesh.data(), flatThisMesh.size(), MPI_FLOAT,
               receiveVerts.data(), receiveVertsSize, MPI_FLOAT,
               0, MPI_COMM_WORLD);

    MPI_Gather(flatThisTris.data(), flatThisTris.size(), my_MPI_SIZE_T,
               receiveTris.data(), receiveTrisSize, my_MPI_SIZE_T,
               0, MPI_COMM_WORLD);

    int gatherSize = numVerts.size();
    std::vector<util::TriangleMesh<float> > meshes;

    // Receive the data and put back into TriangleMesh classes
    if(pid == 0)
    {
        for(int processIdx = 0; processIdx != gatherSize; ++processIdx)
        {
            size_t vIdx = processIdx * receiveVertsSize;
            size_t tIdx = processIdx * receiveTrisSize;

            std::vector<std::array<float, 3> > points(numVerts[processIdx]);
            std::vector<std::array<float, 3> > normals(numVerts[processIdx]);
            std::vector<std::array<size_t, 3> > indexTriangles(numTris[processIdx]);

            for(int i = 0; i != points.size(); ++i)
            {
                points[i][0] = receiveVerts[vIdx++];
                points[i][1] = receiveVerts[vIdx++];
                points[i][2] = receiveVerts[vIdx++];
                normals[i][0] = receiveVerts[vIdx++];
                normals[i][1] = receiveVerts[vIdx++];
                normals[i][2] = receiveVerts[vIdx++];
            }

            for(std::array<size_t, 3>& tri: indexTriangles)
            {
                tri[0] = receiveTris[tIdx++];
                tri[1] = receiveTris[tIdx++];
                tri[2] = receiveTris[tIdx++];
            }

            meshes.emplace_back(points, normals, indexTriangles);
        }
    }

    return meshes;
}

}

#endif
