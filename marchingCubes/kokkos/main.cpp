/*
 * kokkos/main.cpp
 *
 *  Created on: Feb 3, 2017
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

#include <array>
#include <vector>
#include <unordered_map>

#include <string>
#include <string.h>
#include <cstdlib>

#include <ctime>
#include <chrono>
#include <iomanip>

#include "../util/Image3D.h"
#include "../util/TriangleMesh.h"

#include "../util/LoadImage.h"
#include "../util/SaveTriangleMesh.h"

#include "../util/util.h"     // for findCaseId and interpolate
#include "../util/MarchingCubesTables.h"
#include "../util/Errors.h"

#include "../util/Timer.h"
#include "../mantevoCommon/YAML_Doc.hpp"

#include <Kokkos_Core.hpp>

using std::size_t;

template <typename T>
struct MarchingCubesFunctor
{
    MarchingCubesFunctor(
        util::Image3D<T> const& image,
        T const& isoval,
        size_t const& grainDim,
        size_t const& guessSizePoints = 0,
        size_t const& guessSizeTris = 0)
      : image(image), isoval(isoval)
    {
        size_t xBeginIdx = image.xBeginIdx();
        size_t yBeginIdx = image.yBeginIdx();
        size_t zBeginIdx = image.zBeginIdx();

        size_t xEndIdxExtent = image.xEndIdx();
        size_t yEndIdxExtent = image.yEndIdx();
        size_t zEndIdxExtent = image.zEndIdx();

        size_t numSectX = (xEndIdxExtent - xBeginIdx + grainDim - 1) / grainDim;
        size_t numSectY = (yEndIdxExtent - yBeginIdx + grainDim - 1) / grainDim;
        size_t numSectZ = (zEndIdxExtent - zBeginIdx + grainDim - 1) / grainDim;

        numSections = numSectX * numSectY * numSectZ;
        size_t numSectionsPerPage = numSectX * numSectY;

        // TODO derive a neater upper bound. This is just made a guess
        if(guessSizePoints == 0)
        {
            sizePoints = std::max(6*numSections,
                              (xEndIdxExtent - xBeginIdx) *
                              (yEndIdxExtent - yBeginIdx) *
                              (zEndIdxExtent - zBeginIdx) / 10) / numSections;
        }
        else
        {
            sizePoints = (guessSizePoints + numSections - 1) / numSections;
        }

        if(guessSizeTris == 0)
        {
            sizeTris = sizePoints * 2;
        }
        else
        {
            sizeTris = (guessSizeTris + numSections - 1) / numSections;
        }

        positions = Kokkos::View<size_t*[2]>("positions", numSections);
        ptNors = Kokkos::View<T*[6]>("ptNors", numSections * sizePoints);
        tris = Kokkos::View<size_t*[3]>("tris", numSections * sizeTris);
        edgeMap = Kokkos::View<size_t*>("edgeMap", numSections * sizePoints);

        begEndInfo = std::vector<std::array<size_t, 6> >(numSections);
        for(size_t i = 0; i < numSections; ++i)
        {
            size_t xSectIdx = (i % numSectionsPerPage) % numSectX;
            size_t ySectIdx = (i % numSectionsPerPage) / numSectX;
            size_t zSectIdx = (i / numSectionsPerPage);

            size_t xbeg = xSectIdx * grainDim;
            size_t ybeg = ySectIdx * grainDim;
            size_t zbeg = zSectIdx * grainDim;

            size_t xend = std::min(xbeg + grainDim, xEndIdxExtent);
            size_t yend = std::min(ybeg + grainDim, yEndIdxExtent);
            size_t zend = std::min(zbeg + grainDim, zEndIdxExtent);

            positions(i, 0) = i * sizePoints;
            positions(i, 1) = i * sizeTris;

            begEndInfo[i] = {xbeg, ybeg, zbeg, xend, yend, zend};
        }
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator()(int const& sectionId) const // const but this function modifies
                                           // values pointed to by positions,
                                           // ptNors and tris
    {
        // This pointMap will only be used on this sectionId. It will map
        // from global edge indices to indices of ptNors
        std::unordered_map<size_t, size_t> pointMap;

        // Take ptIdx and triIdx by reference. As the algorithm goes,
        // these values will be incremented.
        size_t& ptIdx = positions(sectionId, 0);
        size_t& triIdx = positions(sectionId, 1);

        size_t maxPtIdx = (sectionId + 1) * sizePoints;
        size_t maxTriIdx = (sectionId + 1) * sizeTris;

        pointMap.reserve(maxPtIdx - ptIdx);

        // Retrieve info for this sectionId
        size_t xbeg = begEndInfo[sectionId][0];
        size_t ybeg = begEndInfo[sectionId][1];
        size_t zbeg = begEndInfo[sectionId][2];

        size_t xend = begEndInfo[sectionId][3];
        size_t yend = begEndInfo[sectionId][4];
        size_t zend = begEndInfo[sectionId][5];

        // For each cube, determine whether or not the isosurface insersects
        // the given cube. If so, call processOneCube.
        for(size_t zidx = zbeg; zidx != zend; ++zidx)
        {
            for(size_t yidx = ybeg; yidx != yend; ++yidx)
            {
                auto buffer = image.createBuffer(xbeg, yidx, zidx);
                for(size_t xidx = xbeg; xidx != xend; ++xidx)
                {
                    std::array<T, 8> cubeVertexVals =
                        buffer.getCubeVertexValues(xidx);

                    int cellCaseId = util::findCaseId(cubeVertexVals, isoval);

                    if (cellCaseId == 0 || cellCaseId == 255)
                    {
                        continue;
                    }

                    int newTris = util::numberOfTriangles[cellCaseId];

                    // It could be the case that marchingCubes did not
                    // allocate enough memory. Check that there will
                    // be enough memory. Otherwise, throw an error.
                    if(ptIdx + 3*newTris > maxPtIdx)
                    {
                        std::string e1 = "Try passing a value larger than ";
                        std::string e2 = std::to_string(sizePoints * numSections);
                        std::string e3 = " to the `points_allocate` "
                                         "command line argument";
                        throw util::not_enough_memory_allocated(
                                (e1 + e2 + e3).c_str());
                    }

                    if(triIdx + newTris > maxTriIdx)
                    {
                        std::string e1 = "Try passing a value larger than ";
                        std::string e2 = std::to_string(sizeTris * numSections);
                        std::string e3 = " to the `triangles_allocate` "
                                         "command line argument";
                        throw util::not_enough_memory_allocated(
                                (e1 + e2 + e3).c_str());
                    }

                    // The isosurface intersects this cube.
                    // This function will fill out ptNors, tris and edgeMap
                    // as well as increment ptIdx and triIdx.
                    processOneCube(xidx, yidx, zidx,
                                   sectionId,
                                   cubeVertexVals,
                                   cellCaseId,
                                   pointMap,       // modifies
                                   ptIdx,          // modifies
                                   triIdx);        // modifies
                }
            }
        }
    }

private:
    inline void
    processOneCube(
        size_t const& xidx, size_t const& yidx, size_t const& zidx,
        size_t const&                        sectionId,
        std::array<T, 8> const&              cubeVertexVals,
        int const&                           cellCaseId,
        std::unordered_map<size_t, size_t>&  pointMap, // by non-const reference
        size_t&                              ptIdx,    // by non-const reference
        size_t&                              triIdx    // by non-const reference
        ) const
    {
        // Find the cube configuration of cellCaseId from a lookup table.
        // Then add the triangles of the cube configuration to the
        // ptNors and tris map.

        const int *triEdges = util::caseTrianglesEdges[cellCaseId];

        std::array<std::array<T, 3>, 8> posCube =
            image.getPosCube(xidx, yidx, zidx);
        std::array<std::array<T, 3>, 8> gradCube =
            image.getGradCube(xidx, yidx, zidx);

        for(; *triEdges != -1; triEdges += 3)
        {
            for(int i = 0; i != 3; ++i)
            {
                size_t globalEdgeIndex =
                    image.getGlobalEdgeIndex(xidx, yidx, zidx, triEdges[i]);
                tris(triIdx, i) = globalEdgeIndex;

                if(pointMap.find(globalEdgeIndex) == pointMap.end())
                {
                    pointMap[globalEdgeIndex] = ptIdx;
                    edgeMap[ptIdx] = globalEdgeIndex;

                    const int *vs = util::edgeVertices[triEdges[i]];
                    int v1 = vs[0];
                    int v2 = vs[1];
                    T w = (isoval - cubeVertexVals[v1]) /
                        (cubeVertexVals[v2] - cubeVertexVals[v1]);

                    std::array<T, 3> newPt =
                        util::interpolate(posCube[v1], posCube[v2], w);
                    std::array<T, 3> newNorm =
                        util::interpolate(gradCube[v1], gradCube[v2], w);

                    ptNors(ptIdx, 0) = newPt[0];
                    ptNors(ptIdx, 1) = newPt[1];
                    ptNors(ptIdx, 2) = newPt[2];
                    ptNors(ptIdx, 3) = newNorm[0];
                    ptNors(ptIdx, 4) = newNorm[1];
                    ptNors(ptIdx, 5) = newNorm[2];

                    ++ptIdx;
                }
            }
            ++triIdx;
        }
    }

public:
    // Generate the util::TriangleMesh from the ptNors, tris and
    // edgeMap views.
    util::TriangleMesh<T>
    createTriangleMesh() const
    {
        std::vector<std::array<T, 3> > points;
        std::vector<std::array<T, 3> > normals;
        std::vector<std::array<size_t, 3> > indexTriangles;

        std::unordered_map<size_t, size_t> pointMap;

        // Determine how much space to reserve
        size_t numPts = 0;
        size_t numTris = 0;
        for(size_t sectionId = 0; sectionId != numSections; ++sectionId)
        {
            size_t ptIdx = sectionId * sizePoints;
            size_t endPtIdx = positions(sectionId, 0);
            numPts += (endPtIdx - ptIdx);

            size_t triIdx = sectionId * sizeTris;
            size_t endTriIdx = positions(sectionId, 1);
            numTris += (endTriIdx - triIdx);
        }

        points.reserve(numPts);
        normals.reserve(numPts);
        indexTriangles.reserve(numTris);
        pointMap.reserve(numPts);

        // Add ptNors information to points and normals. Fill out
        // pointMap to remove duplicates. Then fill out the index
        // triangles.
        size_t newPtIdx = 0;
        for(size_t sectionId = 0; sectionId != numSections; ++sectionId)
        {
            size_t ptIdx = sectionId * sizePoints;
            size_t endPtIdx = positions(sectionId, 0);
            for(; ptIdx != endPtIdx; ++ptIdx)
            {
                size_t const& globalEdgeIndex = edgeMap(ptIdx);

                if(pointMap.find(globalEdgeIndex) == pointMap.end())
                {
                    pointMap[globalEdgeIndex] = newPtIdx++;
                    points.push_back(  {ptNors(ptIdx, 0),
                                        ptNors(ptIdx, 1),
                                        ptNors(ptIdx, 2)}  );
                    normals.push_back( {ptNors(ptIdx, 3),
                                        ptNors(ptIdx, 4),
                                        ptNors(ptIdx, 5)}  );
                }
            }

            size_t triIdx = sectionId * sizeTris;
            size_t endTriIdx = positions(sectionId, 1);
            for(; triIdx != endTriIdx; ++triIdx)
            {
                // pointMap maps global edge indices to indices
                // in points and normals.
                indexTriangles.push_back(
                    { pointMap[tris(triIdx, 0)],
                      pointMap[tris(triIdx, 1)],
                      pointMap[tris(triIdx, 2)] } );
            }
        }

        return util::TriangleMesh<T>(points, normals, indexTriangles);
    }

    size_t
    numberOfSections() const
    {
        return numSections;
    }

    double
    pointsPercentage() const
    {
        double percent = 0.0;
        for(size_t sectionId = 0; sectionId != numSections; ++sectionId)
        {
            percent = std::max(
                    percent,
                    (positions(sectionId, 0) - sectionId * sizePoints)
                        / double(sizePoints));
        }
        return percent;
    }

    double
    trianglesPercentage() const
    {
        double percent = 0.0;
        for(size_t sectionId = 0; sectionId != numSections; ++sectionId)
        {
            percent = std::max(
                    percent,
                    (positions(sectionId, 1) - sectionId * sizeTris)
                        / double(sizeTris));
        }
        return percent;
    }

private:
    util::Image3D<T> const& image;
    T const& isoval;

    std::vector<std::array<size_t, 6> > begEndInfo;

    size_t numSections;
    size_t sizePoints;
    size_t sizeTris;

    Kokkos::View<size_t*[2]> positions;
    Kokkos::View<T*[6]> ptNors;
    Kokkos::View<size_t*[3]> tris;
    Kokkos::View<size_t*> edgeMap;
};


int main(int argc, char* argv[])
{
    float isoval;
    bool isovalSet = false;
    char* vtkFile = NULL;
    char* outFile = NULL;
    std::string yamlDirectory = "";
    std::string yamlFileName  = "";

    // numPointsAllocate is the minimum number of points that should be allocated
    // in MarchingCubesFunctor and numTrisAllocate is the minimum number of triangles
    // that should be allocated. If not enough memory is allocated, an error will
    // be thrown. If left to 0, a default value will be used.
    size_t numPointsAllocate = 0;
    size_t numTrisAllocate = 0;

    // To control the granularity of the parallel execution, grainDim is passed
    // to the algorithm. grainDim is the largest number of cubes to be
    // processed in each dimension.
    std::size_t grainDim = 256;

    // Read command line arguments
    for(int i=0; i<argc; i++)
    {
        if( (strcmp(argv[i], "-i") == 0) || (strcmp(argv[i], "-input_file") == 0))
        {
            vtkFile = argv[++i];
        }
        else if( (strcmp(argv[i], "-o") == 0) || (strcmp(argv[i], "-output_file") == 0))
        {
            outFile = argv[++i];
        }
        else if( (strcmp(argv[i], "-v") == 0) || (strcmp(argv[i], "-isoval") == 0))
        {
            isovalSet = true;
            isoval = atof(argv[++i]);
        }
        else if( (strcmp(argv[i], "-g") == 0) || (strcmp(argv[i], "-grain_dim") == 0))
        {
            grainDim = std::stoul(argv[++i]);
        }
        else if( (strcmp(argv[i], "-p") == 0) || (strcmp(argv[i], "-points_allocate") == 0))
        {
            numPointsAllocate = std::stoul(argv[++i]);
        }
        else if( (strcmp(argv[i], "-t") == 0) || (strcmp(argv[i], "-triangles_allocate") == 0))
        {
            numTrisAllocate = std::stoul(argv[++i]);
        }
        else if( (strcmp(argv[i], "-y") == 0) || (strcmp(argv[i], "-yaml_output_file") == 0))
        {
            std::string wholeFile(argv[++i]);

            std::size_t pos = wholeFile.rfind("/");
            if(pos == std::string::npos)
            {
                yamlDirectory = "./";
                yamlFileName = wholeFile;
            }
            else
            {
                yamlDirectory = wholeFile.substr(0, pos + 1);
                yamlFileName = wholeFile.substr(pos + 1);
            }
        }
        else if( (strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0))
        {
            std::cout <<
                "Serial Marching Cubes Options:"  << std::endl <<
                "  -input_file (-i)"              << std::endl <<
                "  -output_file (-o)"             << std::endl <<
                "  -isoval (-v)"                  << std::endl <<
                "  -grain_dim (-g), default 256"  << std::endl <<
                "  -points_allocate (-p)"         << std::endl <<
                "  -triangles_allocate (-t)"      << std::endl <<
                "  -yaml_output_file (-y)"        << std::endl <<
                "  -help (-h)"                    << std::endl;
            return 0;
        }
    }

    if(isovalSet == false || vtkFile == NULL || outFile == NULL)
    {
        std::cout << "Error: isoval, input_file and output_file must be set." << std::endl <<
                     "Try -help" << std::endl;
        return 0;
    }

    // Initialize kokkos
    Kokkos::initialize(argc, argv);

    // Create a yamlDoc. If yamlDirectory and yamlFileName weren't assigned,
    // YAML_Doc will create a file at in the current directory with a
    // timestamp on it.
    YAML_Doc doc("Marching Cubes", "0.2", yamlDirectory, yamlFileName);

    // Add information related to this run to doc.
    doc.add("Marching Cubes Algorithm", "openmp");
    doc.add("Volume image data file path", vtkFile);
    doc.add("Polygonal mesh output file", outFile);
    doc.add("Isoval", isoval);
    doc.add("Grain Dimensions", grainDim);

    // Load the image file
    util::Image3D<float> image = util::loadImage<float>(vtkFile);

    doc.add("File x-dimension", image.xdimension());
    doc.add("File y-dimension", image.ydimension());
    doc.add("File z-dimension", image.zdimension());

    // The number of sections along the x, y and z directions is calculated
    // in MarchingCubes but this is done here so it can be output via doc.
    std::size_t numSectX = (image.xEndIdx() - image.xBeginIdx() + grainDim - 1)
        / grainDim;
    std::size_t numSectY = (image.yEndIdx() - image.xBeginIdx() + grainDim - 1)
        / grainDim;
    std::size_t numSectZ = (image.zEndIdx() - image.zBeginIdx() + grainDim - 1)
        / grainDim;
    doc.add("Number of sections", numSectX * numSectY * numSectZ);

    // Time the output. Timer's constructor starts timing.
    util::Timer runTime;

    // The MarchingCubesFunctor runs the algorithm. The constructor
    // divides the works into sections. As inputs it takes the image
    // loaded at vtkFile and the isoval of the surface to approximate.
    // In addition, the constructor allocates memory in the Kokkos
    // views that will contain the mesh information. numPointsAllocate
    // and numTrisAllocate define how approcimttely how many points to
    // allocate for. If they are zero, marchingCubes makes a guess.
    MarchingCubesFunctor<float> marchingCubes(image, isoval, grainDim,
                                              numPointsAllocate,
                                              numTrisAllocate);

    size_t numSections = marchingCubes.numberOfSections();

    // Use Kokkos and marchingCubes to run the parallel part of the
    // algorithm.
    Kokkos::parallel_for(numSections, marchingCubes);

    // Kokkos may return before the for loop is complete.
    // Wait for all the threads to be done here.
    Kokkos::fence();

    // Now that the data has been collected in parallel, merge it and output
    // the final mesh. This process ensures that no duplicate points are
    // output.
    util::TriangleMesh<float> polygonalMesh = marchingCubes.createTriangleMesh();

    // End timing
    runTime.stop();

    // Report mesh information
    doc.add("Number of vertices in mesh", polygonalMesh.numberOfVertices());
    doc.add("Number of triangles in mesh", polygonalMesh.numberOfTriangles());

    // Report the max percentage of memory allocation used by the threads
    doc.add("Max % memory allocation for points",
            marchingCubes.pointsPercentage());
    doc.add("Max % memory allocation for triangles",
            marchingCubes.trianglesPercentage());

    // Report timing information
    doc.add("Total Program CPU Time (clicks)", runTime.getTotalTicks());
    doc.add("Total Program CPU Time (seconds)", runTime.getCPUtime());
    doc.add("Total Program WALL Time (seconds)", runTime.getWallTime());

    // Generate the YAML file. The file will be both saved and printed to console.
    std::cout << doc.generateYAML();

    // Save the polygonal mesh to the output file
    util::saveTriangleMesh(polygonalMesh, outFile);

    // Don't forget to tell Kokkos bye!
    Kokkos::finalize();
}
