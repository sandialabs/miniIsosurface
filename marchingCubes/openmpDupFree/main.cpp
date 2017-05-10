/*
 * openmpDupFree/main.cpp
 *
 *  Created on: Jan 25, 2017
 *      Author: dbourge
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

#include <omp.h>

#include "../util/Image3D.h"
#include "../util/TriangleMesh.h"

#include "../util/LoadImage.h"
#include "../util/SaveTriangleMesh.h"

#include "../util/util.h"     // for findCaseId and interpolate
#include "../util/MarchingCubesTables.h"
#include "../util/DuplicateRemover.h"

#include "../util/Timer.h"
#include "../mantevoCommon/YAML_Doc.hpp"

using std::size_t;

template <typename T>
void
sectionOfMarchingCubes(
    size_t const& xbeg, size_t const& ybeg, size_t const& zbeg,
    size_t const& xend, size_t const& yend, size_t const& zend,
    T const&                                isoval,
    util::Image3D<T> const&                 image,
    std::vector<std::array<T, 3> >&         points,         // reference
    std::vector<std::array<T, 3> >&         normals,        // reference
    std::vector<std::array<size_t, 3> >&  indexTriangles, // reference
    std::unordered_map<size_t, size_t>& pointMap)       // reference
{
    // For each cube, determine whether or not the isosurface intersects
    // the given cube. If so, first find the cube configuration from a lookup
    // table. Then add the triangles of that cube configuration to points,
    // normals and triangles. Use pointMap to not add any duplicates on
    // this thread to points or normals

    size_t ptIdx = points.size();
    for(size_t zidx = zbeg; zidx != zend; ++zidx)
    {
        for(size_t yidx = ybeg; yidx != yend; ++yidx)
        {
            // A buffer is used to improve cache efficency when retrieving
            // vertex values of each cube.
            auto buffer = image.createBuffer(xbeg, yidx, zidx);

            for(size_t xidx = xbeg; xidx != xend; ++xidx)
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
                    std::array<size_t, 3> tri;
                    for(int i = 0; i != 3; ++i)
                    {
                        size_t globalEdgeIndex =
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
util::TriangleMesh<T>
MarchingCubes(util::Image3D<T> const& image, T const& isoval,
    size_t const& nSectionsX, size_t const& nSectionsY, size_t const& nSectionsZ)
{
    // The marching cubes algorithm creates a polygonal mesh to approximate an
    // isosurface from a three-dimensional discrete scalar field.

    // The three-dimensional discrete scale field is stored in image and the
    // constant value of the isosurface to approximate is isoval.

    // The polygonal mesh is represented by the TriangleMesh class.
    // A vector of index triangles and the corresponding points and normal
    // vectors of each triangle's vertices are needed to create an instance
    // of TriangleMesh.
    std::vector<std::array<T, 3> > points;
    std::vector<std::array<T, 3> > normals;
    // Example:
    //   If indexTriangles[5] = {1, 4, 2}, then the polygonal mesh will contain
    //   a triangle with vertices at points[1], points[4] and points[2].
    std::vector<std::array<size_t, 3> > indexTriangles;

    // If multiple threads or processes are present, points may be added
    // more than once. To fix this, duplicateTracker will be filled with pairs
    // containing the point index in points/normals and the global edge index.
    // Later, util::duplicateRemover will remove the duplicates and return
    // the final, duplicate free mesh.
    std::vector<std::pair<size_t, size_t> > duplicateTracker;

    // Using OpenMP, this code is ran in parallel one section at a time.
    // The size of each section is determined by grainDim, which is the
    // largest number of cubes in each dimenstion.
    // Example:
    //   suppose image is 501x501x501 (so the algorithm will look
    //   at 500x500x500 cubes)
    //   suppose grainDim is 300
    //   Then there will be 8 partitions:
    //     xs, ys, zs
    //     0-299, 0-299, 0-299
    //     0-299, 0-299, 300-499
    //     0-299, 300-499, 300-499
    //     ...

    size_t xBeginIdx = image.xBeginIdx();
    size_t yBeginIdx = image.yBeginIdx();
    size_t zBeginIdx = image.zBeginIdx();

    size_t xEndIdxExtent = image.xEndIdx();
    size_t yEndIdxExtent = image.yEndIdx();
    size_t zEndIdxExtent = image.zEndIdx();

    // The util::Indexer statically distributes the elements that the
    // section is in charge of.
    // The operator function of a util::Indexer takes a section index
    // and returns the corresponding index starting at that section.
    // Each section W has total_elements_on_section / n_sections or
    // total_elements_on_section / n_sections + 1 elements in it such
    // that indexer(n_sections) = total_elements_on_section.
    util::Indexer indexerX(xEndIdxExtent - xBeginIdx, nSectionsX);
    util::Indexer indexerY(yEndIdxExtent - yBeginIdx, nSectionsY);
    util::Indexer indexerZ(zEndIdxExtent - zBeginIdx, nSectionsZ);

    size_t nSections = nSectionsX * nSectionsY * nSectionsZ;
    size_t nSectionsPerPage = nSectionsX * nSectionsY;

    #pragma omp parallel
    {
        // Each openMP thread manages it's own threadPoints, threadNormals,
        // threadIndexTriangles and threadPointMap.
        std::vector<std::array<T, 3> > threadPoints;
        std::vector<std::array<T, 3> > threadNormals;
        std::vector<std::array<size_t, 3> > threadIndexTriangles;

        // threadPointMap will be a dictionary from global edge indices
        // to indices in threadPoints and threadNormals. Note that
        // threadPoints and threadNormals share the same indices
        // but are only unique among this section/thread of execution.
        std::unordered_map<size_t, size_t> threadPointMap;

        #pragma omp for nowait
        for(size_t i = 0; i < nSections; ++i)
        {
            // Determine the coordinates of this section.
            size_t xSectIdx = (i % nSectionsPerPage) % nSectionsX;
            size_t ySectIdx = (i % nSectionsPerPage) / nSectionsX;
            size_t zSectIdx = (i / nSectionsPerPage);

            size_t xbeg = indexerX(xSectIdx);
            size_t ybeg = indexerY(ySectIdx);
            size_t zbeg = indexerZ(zSectIdx);

            // indexerW(nSectionsW) == wEndIdxExtent - wBeginIdx
            size_t xend = indexerX(xSectIdx + 1);
            size_t yend = indexerY(ySectIdx + 1);
            size_t zend = indexerZ(zSectIdx + 1);

            // How does this work? TODO
            // For performance reasons, rehashing the unordered map.
            size_t approxNumberOfEdges = 3*(xend-xbeg)*(yend-ybeg)*(zend-zbeg);
            size_t mapSize = approxNumberOfEdges / 8 + 6;
            threadPointMap.rehash(mapSize);

            // Variables for this thread of execution are given by reference and
            // will be modified.
            sectionOfMarchingCubes(
                xbeg, ybeg, zbeg,           // constant inputs
                xend, yend, zend,           // constant inputs
                isoval, image,              // constant inputs
                threadPoints,              // for modification, taken by reference
                threadNormals,             // for modification, taken by reference
                threadIndexTriangles,      // for modification, taken by reference
                threadPointMap);           // for modification, taken by reference
        }

        // As each section is complete, the mesh information is added to
        // points, normals and indexTriangles and duplicateTracker. The triangles
        // in threadIndexTriangles will have to be offset to have indices
        // that coincide with points and not threadPoints.
        //
        // critical ensures that this block of code will execute one
        // thread at a time.
        #pragma omp critical
        {
            // Reserve space in points, normals and indexTriangles so that
            // these vectors don't have to resize themselves constantly.
            points.reserve(points.size() + threadPoints.size());
            normals.reserve(normals.size() + threadNormals.size());
            indexTriangles.reserve(
                indexTriangles.size() + threadIndexTriangles.size());

            size_t offset = points.size();

            points.insert(
                points.end(), threadPoints.begin(), threadPoints.end());
            normals.insert(
                normals.end(), threadNormals.begin(), threadNormals.end());

            // The points refered to in each tri need to refer to the points
            // in the points vector, not threadPoints.
            for(std::array<size_t, 3>& tri: threadIndexTriangles)
            {
                tri[0] = tri[0] + offset;
                tri[1] = tri[1] + offset;
                tri[2] = tri[2] + offset;
            }

            indexTriangles.insert(
                indexTriangles.end(),
                threadIndexTriangles.begin(), threadIndexTriangles.end());

            // duplicateTracker provides information for the util::duplicateRemover
            // function.
            duplicateTracker.resize(offset + threadPointMap.size());
            for(auto iter = threadPointMap.begin(); iter != threadPointMap.end(); ++iter)
            {
                size_t const& pointIndex = offset + iter->second;
                size_t const& globalEdgeIndex = iter->first;

                duplicateTracker[pointIndex] =
                    std::make_pair(pointIndex, globalEdgeIndex);
            }
        }
    }

    // If there are no triangles whatsover, return here.
    if(indexTriangles.size() == 0)
    {
        return util::TriangleMesh<T>();
    }

    return util::duplicateRemover(
        duplicateTracker, points, normals, indexTriangles);
}

int main(int argc, char* argv[])
{
    float isoval;
    bool isovalSet = false;
    char* vtkFile = NULL;
    char* outFile = NULL;
    std::string yamlDirectory = "";
    std::string yamlFileName  = "";

    // To control the granularity of the parallel execution,
    // specify how many sections should be in the X, Y and Z direction.
    std::size_t nSectionsX = 1;
    std::size_t nSectionsY = 1;

    std::size_t nSectionsZ;
    // set nSectionsZ to the number of openmp threads being run
    // Calling omp_get_num_threads() outside of a parallel directive
    // will return 1
    #pragma omp parallel
    {
        if(omp_get_thread_num() == 0)
        {
            nSectionsZ = omp_get_num_threads();
        }
    }

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
        else if( (strcmp(argv[i], "-sx") == 0) || (strcmp(argv[i], "-sections_x") == 0))
        {
            nSectionsX = std::stoul(argv[++i]);
        }
        else if( (strcmp(argv[i], "-sy") == 0) || (strcmp(argv[i], "-sections_y") == 0))
        {
            nSectionsY = std::stoul(argv[++i]);
        }
        else if( (strcmp(argv[i], "-sz") == 0) || (strcmp(argv[i], "-sections_z") == 0))
        {
            nSectionsZ = std::stoul(argv[++i]);
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
                "  -sections_x (-sx)"             << std::endl <<
                "  -sections_y (-sy)"             << std::endl <<
                "  -sections_z (-sz)"             << std::endl <<
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

    // Create a yamlDoc. If yamlDirectory and yamlFileName weren't assigned,
    // YAML_Doc will create a file at in the current directory with a
    // timestamp on it.
    YAML_Doc doc("Marching Cubes", "0.2", yamlDirectory, yamlFileName);

    // Add information related to this run to doc.
    doc.add("Marching Cubes Algorithm", "openmpDupFree");
    doc.add("Volume image data file path", vtkFile);
    doc.add("Polygonal mesh output file", outFile);
    doc.add("Isoval", isoval);

    // Load the image file
    util::Image3D<float> image = util::loadImage<float>(vtkFile);

    // Readjust nSections if they are set too large.
    if(nSectionsX > image.xdimension() - 1)
        nSectionsX = image.xdimension() - 1;
    if(nSectionsY > image.ydimension() - 1)
        nSectionsY = image.ydimension() - 1;
    if(nSectionsZ > image.zdimension() - 1)
        nSectionsZ = image.zdimension() - 1;

    doc.add("Number of X sections", nSectionsX);
    doc.add("Number of Y sections", nSectionsY);
    doc.add("Number of Z sections", nSectionsZ);
    doc.add("Number of sections", nSectionsX * nSectionsY * nSectionsZ);

    doc.add("File x-dimension", image.xdimension());
    doc.add("File y-dimension", image.ydimension());
    doc.add("File z-dimension", image.zdimension());

    // Time the output. Timer's constructor starts timing.
    util::Timer runTime;

    // MarchingCubes runs the algorithm. As inputs it takes the image
    // loaded at vtkFile and the isoval of the surface to approximate. It's
    // output is a TriangleMesh which stores the mesh as a vector
    // of triangles.
    util::TriangleMesh<float> polygonalMesh =
        MarchingCubes(image, isoval, nSectionsX, nSectionsY, nSectionsZ);

    // End timing
    runTime.stop();

    // Report mesh information
    doc.add("Number of vertices in mesh", polygonalMesh.numberOfVertices());
    doc.add("Number of triangles in mesh", polygonalMesh.numberOfTriangles());

    // Report timing information
    doc.add("Total Program CPU Time (clicks)", runTime.getTotalTicks());
    doc.add("Total Program CPU Time (seconds)", runTime.getCPUtime());
    doc.add("Total Program WALL Time (seconds)", runTime.getWallTime());

    // Generate the YAML file. The file will be both saved and printed to console.
    std::cout << doc.generateYAML();

    // Save the polygonal mesh to the output file
    util::saveTriangleMesh(polygonalMesh, outFile);
}
