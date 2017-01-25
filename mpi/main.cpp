/*
 * mpi/main.cpp
 *
 *  Created on: Jan 25, 2017
 *      Author: dbourge
 */

#include <array>
#include <vector>
#include <unordered_map>

#include <cstdlib>

#include <ctime>
#include <chrono>
#include <iomanip>

#include "../util/util.h"     // for findCaseId and interpolate
#include "../util/MarchingCubesTables.h"
#include "../util/TypeInfo.h"

#include "../util/LoadImage.h"
#include "../util/SaveTriangleMesh.h"

#include "../util/Image3D.h"
#include "../util/TriangleMesh.h"

#include "mpi.h"

template <typename T>
void
sectionOfMarchingCubes(
    T const&                                isoval,
    util::Image3D<T> const&                 image,
    std::vector<std::array<T, 3> >&         points,         // reference
    std::vector<std::array<T, 3> >&         normals,        // reference
    std::vector<std::array<unsigned, 3> >&  indexTriangles, // reference
    std::unordered_map<unsigned, unsigned>& pointMap)       // reference
{
    // For each cube, determine whether or not the isosurface intersects
    // the given cube. If so, first find the cube configuration from a lookup
    // table. Then add the triangles of that cube configuration to points,
    // normals and triangles. Use pointMap to not add any duplicate to
    // points or normals.

    unsigned xbeg = image.xBeginIdx();
    unsigned ybeg = image.yBeginIdx();
    unsigned zbeg = image.zBeginIdx();

    unsigned xend = image.xEndIdx();
    unsigned yend = image.yEndIdx();
    unsigned zend = image.zEndIdx();

    unsigned ptIdx = points.size();
    for(unsigned zidx = zbeg; zidx != zend; ++zidx)
    {
        for(unsigned yidx = ybeg; yidx != yend; ++yidx)
        {
            // A buffer is used to improve cache efficency when retrieving
            // vertex values of each cube.
            auto buffer = image.createBuffer(xbeg, yidx, zidx);

            for(unsigned xidx = xbeg; xidx != xend; ++xidx)
            {
                // For each x, y, z index, get the corresponding 8 scalar values
                // of a cube with one corner being at the x, y, z index.
                std::array<T, 8> cubeVertexVals = buffer.getCubeVertexValues(xidx);

                int cellCaseId = util::findCaseId(cubeVertexVals, isoval);

                // The 0th and 255th cellCaseId are the cases where the isosurface
                // does not intersect the cube.
                if (cellCaseId == 0 || cellCaseId == 255)
                {
                    continue;
                }

                // From a pre-generated list of possible cube configurations,
                // get the corresponding "triangles" to this caseId. Each
                // triangle contains 3 cube edge endices, where each vertex
                // of the triangle lies uniquely on one of the edges.
                const int *triEdges = util::caseTrianglesEdges[cellCaseId];
                // triEdges has the form [a1, b1, c1, ..., an, bn, cn, -1]
                // where ax, bx, cx are the edge indices of a triangle for this
                // cube configuration. n can be from 1 to 5.

                // No buffer is used for retrieving scalar values and
                // calculating gradients.
                std::array<std::array<T, 3>, 8> posCube =
                    image.getPosCube(xidx, yidx, zidx);
                std::array<std::array<T, 3>, 8> gradCube =
                    image.getGradCube(xidx, yidx, zidx);

                // For each edge in each triangle, find it's global edge index.
                // If the point and normal for that  global edge index has not
                // been calculated, calculate it and add it to the points and
                // normals vector. Plus, keep track of the indices in points
                // and normals that will make up the next index triangle to
                // add to indexTriangles.
                for(; *triEdges != -1; triEdges += 3)
                {
                    // tri contains indices to points and normals
                    std::array<unsigned, 3> tri;
                    for(int i = 0; i != 3; ++i)
                    {
                        unsigned globalEdgeIndex =
                            image.getGlobalEdgeIndex(xidx, yidx, zidx, triEdges[i]);

                        if(pointMap.find(globalEdgeIndex) != pointMap.end())
                        {
                            // This global edge index has an associated index
                            // corresponding to its point and normal in points
                            // and normals.

                            tri[i] = pointMap[globalEdgeIndex];
                        }
                        else
                        {
                            // This point and normal value on this global edge
                            // index has not been calculated.

                            pointMap[globalEdgeIndex] = ptIdx;
                            tri[i] = ptIdx++;

                            const int *vs = util::edgeVertices[triEdges[i]];
                            int v1 = vs[0];
                            int v2 = vs[1];
                            T w = (isoval - cubeVertexVals[v1]) /
                                (cubeVertexVals[v2] - cubeVertexVals[v1]);

                            std::array<T, 3> newPt =
                                util::interpolate(posCube[v1], posCube[v2], w);
                            std::array<T, 3> newNorm =
                                util::interpolate(gradCube[v1], gradCube[v2], w);

                            points.push_back(newPt);
                            normals.push_back(newNorm);
                        }
                    }

                    // TODO Is this check needed?
                    if(tri[0] != tri[1] && tri[1] != tri[2] && tri[2] != tri[0])
                    {
                        indexTriangles.push_back(tri);
                    }
                }
            }
        }
    }
}

template <typename T>
std::vector<T>
readSectionData(
    unsigned xbeg, unsigned ybeg, unsigned zbeg,
    unsigned xend, unsigned yend, unsigned zend,
    const char*                         file,
    util::TypeInfo const&               ti,
    std::array<unsigned, 3> const&      globalDim)
{
    std::ifstream stream(file);
    if (!stream)
        throw util::file_not_found(file);

    // stream is taken by reference. When finished, this function puts stream
    // ahead of all of the header information.
    util::skipHeader(stream);

    unsigned nPointsIgnore =
        xbeg + (ybeg * globalDim[0]) + (zbeg * globalDim[0] * globalDim[1]);

    // stream is taken by reference. Place stream nPointsIgnore ahead where
    // each point is of size ti.size().
    util::streamIgnore(stream, nPointsIgnore, ti.size());

    unsigned nXpoints = xend - xbeg;
    unsigned nYpoints = yend - ybeg;
    unsigned nZpoints = zend - zbeg;
    unsigned nPointsInSection = nXpoints * nYpoints * nZpoints;

    unsigned nXpointsIgnore = globalDim[0] - nXpoints;
    unsigned nYpointsIgnore = globalDim[0] * (globalDim[1] - nYpoints);

    std::size_t readXlineSize = nXpoints * ti.size();
    std::size_t totalReadSize = nPointsInSection * ti.size();

    std::vector<char> rbufRead(totalReadSize);

    unsigned imageDataIdx = 0;
    for(unsigned iZline = 0; iZline < nZpoints; ++iZline)
    {
        for(unsigned iYline = 0; iYline < nYpoints; ++iYline)
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
loadImageSections(const char* file, unsigned const& grainDim)
{
    std::ifstream stream(file);
    if (!stream)
        throw util::file_not_found(file);

    std::array<unsigned, 3> dim;
    std::array<T, 3> spacing;
    std::array<T, 3> zeroPos;
    unsigned npoints;
    util::TypeInfo ti;

    // These variables are all taken by reference
    loadHeader(stream, dim, spacing, zeroPos, npoints, ti);
    stream.close();

    unsigned xBeginIdx = 0;
    unsigned yBeginIdx = 0;
    unsigned zBeginIdx = 0;

    unsigned xEndIdxExtent = dim[0] - 1;
    unsigned yEndIdxExtent = dim[1] - 1;
    unsigned zEndIdxExtent = dim[2] - 1;

    unsigned numSectX = (xEndIdxExtent - xBeginIdx + grainDim - 1) / grainDim;
    unsigned numSectY = (yEndIdxExtent - yBeginIdx + grainDim - 1) / grainDim;
    unsigned numSectZ = (zEndIdxExtent - zBeginIdx + grainDim - 1) / grainDim;

    unsigned numSections = numSectX * numSectY * numSectZ;
    unsigned numSectionsPerPage = numSectX * numSectY;

    int pid = MPI::COMM_WORLD.Get_rank();
    int nProcesses = MPI::COMM_WORLD.Get_size();

    unsigned sectPerProcess = (numSections + nProcesses - 1) / nProcesses;
    unsigned split = nProcesses + numSections - sectPerProcess * nProcesses;

    unsigned startSectNum, endSectNum;
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
    for(unsigned i = startSectNum; i != endSectNum; ++i)
    {
        // Determine the coordinates of this section.
        unsigned xSectIdx = (i % numSectionsPerPage) % numSectX;
        unsigned ySectIdx = (i % numSectionsPerPage) / numSectX;
        unsigned zSectIdx = (i / numSectionsPerPage);

        // Let w be either x, y, or z. wBegIdx/wEndIdx and wDataBeg/wDataEnd
        // will most likely refer to different ranges. wBegIdx/wEndIdx defines
        // the range of indicies that image.createBuffer and image.getGradCube
        // have valid inputs.
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
        unsigned xBegIdx = xSectIdx * grainDim;
        unsigned yBegIdx = ySectIdx * grainDim;
        unsigned zBegIdx = zSectIdx * grainDim;

        unsigned xEndIdx = std::min(xBegIdx + grainDim, xEndIdxExtent);
        unsigned yEndIdx = std::min(yBegIdx + grainDim, yEndIdxExtent);
        unsigned zEndIdx = std::min(zBegIdx + grainDim, zEndIdxExtent);

        unsigned xDataBeg = xBegIdx == 0   ?   0   :   xBegIdx - 1;
        unsigned yDataBeg = yBegIdx == 0   ?   0   :   yBegIdx - 1;
        unsigned zDataBeg = zBegIdx == 0   ?   0   :   zBegIdx - 1;

        unsigned xDataEnd = std::min(xEndIdx + 2, dim[0]);
        unsigned yDataEnd = std::min(yEndIdx + 2, dim[1]);
        unsigned zDataEnd = std::min(zEndIdx + 2, dim[2]);

        std::vector<T> imageData = readSectionData<T>(
            xDataBeg, yDataBeg, zDataBeg, xDataEnd, yDataEnd, zDataEnd,
            file, ti, dim);

        images.emplace_back(
            imageData,
            spacing,
            zeroPos,
            std::array<unsigned, 3>({xBegIdx, yBegIdx, zBegIdx}),
            std::array<unsigned, 3>({xEndIdx, yEndIdx, zEndIdx}),
            std::array<unsigned, 3>({xDataBeg, yDataBeg, zDataBeg}),
            std::array<unsigned, 3>({xDataEnd, yDataEnd, zDataEnd}),
            dim);
    }

    return images;
}

template <typename T>
util::TriangleMesh<T>
MarchingCubes(std::vector<util::Image3D<T> > const& images, T const& isoval)
{
    std::vector<std::array<T, 3> > processPoints;
    std::vector<std::array<T, 3> > processNormals;
    std::vector<std::array<unsigned, 3> > processIndexTriangles;

    std::unordered_map<unsigned, unsigned> processPointMap;

    for(util::Image3D<T> const& image: images)
    {
        sectionOfMarchingCubes(
            isoval, image,                  // constant inputs
            processPoints,                  // for modification, taken by reference
            processNormals,                 // for modification, taken by reference
            processIndexTriangles,          // for modification, taken by reference
            processPointMap);               // for modification, taken by reference
    }

    // TODO The mesh information is not communicated across processes in the original
    // implementation. Instead, each process outputs a process unique mesh which
    // then gets output to n output files, where n is the number of processes.
    //
    // Also, the original implementation tried to remove duplicate points. But because
    // this function is done sequentially on this process and all the point information
    // will be stored in this set of process vectors, there won't be any duplicated
    // points--on this process. However, there could be duplicated points across all of
    // the processes.
    //
    // The TODO is:
    //   -Should duplicate points be removed and there only be one output?
    //   -Or should there be n output files? This is what is implemented now.  In this
    //   case the duplicated points are needed.

    return util::TriangleMesh<T>(processPoints, processNormals, processIndexTriangles);
}

int main(int argc, char* argv[])
{
    MPI::Init(argc, argv);

    char* vtkFile = argv[1];
    std::string outFile = argv[2];
    float isoval = atof(argv[3]);

    // To control the granularity of the parallel execution, grainDim is passed
    // to the algorithm. grainDim is the largest number of cubes to be
    // processed in each dimension.
    unsigned grainDim = 256;
    if(argc == 5)
    {
        grainDim = atoi(argv[4]);
    }

    std::vector<util::Image3D<float> > images =
        loadImageSections<float>(vtkFile, grainDim);

    // Time the output TODO
    // std::clock_t c_start = std::clock();
    // auto t_start = std::chrono::high_resolution_clock::now();

    util::TriangleMesh<float> polygonalMesh = MarchingCubes(images, isoval);

    // End timing
    // std::clock_t c_end = std::clock();
    // auto t_end = std::chrono::high_resolution_clock::now();

    // Print the output
    //std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
    //          << (c_end-c_start) / (1.0 * CLOCKS_PER_SEC) << " s\n"
    //          << "Wall clock time passed: "
    //          << std::chrono::duration<double>(t_end-t_start).count()
    //          << " s\n";

    int pid = MPI::COMM_WORLD.Get_rank();
	outFile = outFile + "." + std::to_string(static_cast<long long int>(pid));
    util::saveTriangleMesh(polygonalMesh, outFile.c_str());

    MPI::Finalize();
}
