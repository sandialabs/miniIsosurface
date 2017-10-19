/*
 * mpi/mpiutil.h
 *
 *  Created on: Feb 2, 2017
 *      Author: dbourge
 *
 * miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
 * See LICENSE.txt for details.
 *
 * Copyright (c) 2017
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
 * the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
 * certain rights in this software.
 */

#ifndef MPIUTHELPERFUNCTIONSM
#define MPIUTHELPERFUNCTIONSM

#include <array>
#include <vector>
#include <unordered_map>

#include <algorithm>

#include "../util/TriangleMesh.h"
#include "../util/util.h"

#include <mpi.h>

namespace mpiutil {

template <typename T>
std::vector<T>
readSectionData(
    size_t xbeg, size_t ybeg, size_t zbeg,
    size_t xend, size_t yend, size_t zend,
    const char*                         file,
    util::TypeInfo const&               ti,
    std::array<size_t, 3> const&        globalDim)
{
    std::ifstream stream(file);
    if (!stream)
        throw util::file_not_found(file);

    // stream is taken by reference. When finished, this function puts stream
    // ahead of all of the header information.
    util::skipHeader(stream);

    size_t nPointsIgnore =
        xbeg + (ybeg * globalDim[0]) + (zbeg * globalDim[0] * globalDim[1]);

    // stream is taken by reference. Place stream nPointsIgnore ahead where
    // each point is of size ti.size().
    util::streamIgnore(stream, nPointsIgnore, ti.size());

    size_t nXpoints = xend - xbeg;
    size_t nYpoints = yend - ybeg;
    size_t nZpoints = zend - zbeg;
    size_t nPointsInSection = nXpoints * nYpoints * nZpoints;

    size_t nXpointsIgnore = globalDim[0] - nXpoints;
    size_t nYpointsIgnore = globalDim[0] * (globalDim[1] - nYpoints);

    std::size_t readXlineSize = nXpoints * ti.size();
    std::size_t totalReadSize = nPointsInSection * ti.size();

    std::vector<char> rbufRead(totalReadSize);

    size_t imageDataIdx = 0;
    for(size_t iZline = 0; iZline < nZpoints; ++iZline)
    {
        for(size_t iYline = 0; iYline < nYpoints; ++iYline)
        {
            stream.read(&rbufRead[imageDataIdx], readXlineSize);
            util::streamIgnore(stream, nXpointsIgnore, ti.size());
            imageDataIdx += readXlineSize;
        }
        util::streamIgnore(stream, nYpointsIgnore, ti.size());
    }

    std::vector<T> data(nPointsInSection);
    util::convertBufferWithTypeInfo(rbufRead.data(), ti, nPointsInSection, data.data());

    return data;
}

template <typename T>
std::vector<util::Image3D<T> >
loadImageSections(const char* file,
    size_t const& nSectionsX, size_t const& nSectionsY, size_t const& nSectionsZ)
{
    std::ifstream stream(file);
    if (!stream)
        throw util::file_not_found(file);

    std::array<size_t, 3> dim;
    std::array<T, 3> spacing;
    std::array<T, 3> zeroPos;
    size_t npoints;
    util::TypeInfo ti;

    // These variables are all taken by reference
    loadHeader(stream, dim, spacing, zeroPos, npoints, ti);
    stream.close();

    size_t xBeginIdx = 0;
    size_t yBeginIdx = 0;
    size_t zBeginIdx = 0;

    size_t xEndIdxExtent = dim[0] - 1;
    size_t yEndIdxExtent = dim[1] - 1;
    size_t zEndIdxExtent = dim[2] - 1;

    util::Indexer indexerX(xEndIdxExtent - xBeginIdx, nSectionsX);
    util::Indexer indexerY(yEndIdxExtent - yBeginIdx, nSectionsY);
    util::Indexer indexerZ(zEndIdxExtent - zBeginIdx, nSectionsZ);

    size_t nSections = nSectionsX * nSectionsY * nSectionsZ;
    size_t nSectionsPerPage = nSectionsX * nSectionsY;

    int pid = MPI::COMM_WORLD.Get_rank();
    int nProcesses = MPI::COMM_WORLD.Get_size();

    size_t sectPerProcess = (nSections + nProcesses - 1) / nProcesses;
    size_t split = nProcesses + nSections - sectPerProcess * nProcesses;

    size_t startSectNum, endSectNum;
    if (pid < split)
    {
        startSectNum = pid * sectPerProcess;
        endSectNum = startSectNum + sectPerProcess;
    }
    else
    {
        startSectNum = split * sectPerProcess +
                       (pid - split) * (sectPerProcess - 1);
        endSectNum = startSectNum + sectPerProcess - 1;
    }

    std::vector<util::Image3D<T> > images;
    for(size_t i = startSectNum; i != endSectNum; ++i)
    {
        // Determine the coordinates of this section.
        size_t xSectIdx = (i % nSectionsPerPage) % nSectionsX;
        size_t ySectIdx = (i % nSectionsPerPage) / nSectionsX;
        size_t zSectIdx = (i / nSectionsPerPage);

        // Let w be either x, y, or z. wBegIdx/wEndIdx and wDataBeg/wDataEnd
        // will most likely refer to different ranges. wBegIdx/wEndIdx defines
        // the range of indicies that image.createBuffer and image.getGradCube
        // have valid inputs, whereas wDataBeg/wDataEnd are ranges of the
        // ghost cells.
        //
        // Calling getGradCube with (xBegIdx, yBegIdx, zEndIdx-1) will try to
        // access vertex value data at (xBegIdx-1, yBegIdx, zEndIdx-1) as well
        // as at (xBegIdx, yBegIdx, zEndIdx+1) even though these values are
        // outside of the [wBegIdx, wEndIdx) ranges. The actual range for data
        // values is instead given by wDataBeg/wDataEnd.
        //
        // This is done so that for the same inputs, the mpi implementaiton will
        // output equivalent meshes as the reference or openmp implementations.
        //
        // Example:
        //  Suppose we have a 1D image that has 100 points and we want to have
        //  1 section. Then for the first section:
        //    xDataBeg = 0, xDataEnd = 100
        //    xBegIdx = 0, xEndIdx = 99
        // Example:
        //  Suppose we have a 1D image that has 512 points and we want to have
        //  2 sections of (close to) the same size
        //   The first section:
        //     xDataBeg = 0, xDataEnd = 258
        //     xBegIdx = 0,  xEndIdx = 256
        //   The second section:
        //     xDataBeg = 255, xDataEnd = 512
        //     xBegIdx = 256, xEndIdx = 511
        // Example:
        //   Suppose we have a 1D image that has 300 points and we want to have
        //   3 sections of (close to) the same size.
        //   The first section:
        //     xDataBeg = 0, xDataEnd = 102
        //     xBegIdx = 0, xEndIdx = 100
        //   The second section:
        //     xDataBeg = 99 , xDataEnd = 202
        //     xBegIdx = 100, xEndIdx = 200
        //   The third section:
        //     xDataBeg = 199, xDataEnd = 300
        //     xBegIdx = 200, xEndIdx = 299
        size_t xBegIdx = indexerX(xSectIdx);
        size_t yBegIdx = indexerY(ySectIdx);
        size_t zBegIdx = indexerZ(zSectIdx);

        size_t xEndIdx = indexerX(xSectIdx + 1);
        size_t yEndIdx = indexerY(ySectIdx + 1);
        size_t zEndIdx = indexerZ(zSectIdx + 1);

        size_t xDataBeg = xBegIdx == 0   ?   0   :   xBegIdx - 1;
        size_t yDataBeg = yBegIdx == 0   ?   0   :   yBegIdx - 1;
        size_t zDataBeg = zBegIdx == 0   ?   0   :   zBegIdx - 1;

        size_t xDataEnd = std::min(xEndIdx + 2, dim[0]);
        size_t yDataEnd = std::min(yEndIdx + 2, dim[1]);
        size_t zDataEnd = std::min(zEndIdx + 2, dim[2]);

        std::vector<T> imageData = readSectionData<T>(
            xDataBeg, yDataBeg, zDataBeg, xDataEnd, yDataEnd, zDataEnd,
            file, ti, dim);

        images.emplace_back(
            imageData,
            spacing,
            zeroPos,
            std::array<size_t, 3>({xBegIdx, yBegIdx, zBegIdx}),
            std::array<size_t, 3>({xEndIdx, yEndIdx, zEndIdx}),
            std::array<size_t, 3>({xDataBeg, yDataBeg, zDataBeg}),
            std::array<size_t, 3>({xDataEnd, yDataEnd, zDataEnd}),
            dim);
    }

    return images;
}


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
