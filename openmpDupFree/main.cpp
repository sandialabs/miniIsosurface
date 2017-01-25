#include <array>
#include <vector>
#include <unordered_map>

#include <algorithm>

#include <cstdlib>

#include <ctime>
#include <chrono>
#include <iomanip>

#include "../util/util.h"     // for findCaseId and interpolate
#include "../util/MarchingCubesTables.h"
#include "../util/DuplicateRemover.h"

#include "../util/LoadImage.h"
#include "../reference/SaveTriangleMesh.h" // Same as at reference TODO

#include "../util/Image3D.h"
#include "../reference/TriangleMesh.h"     // Same as at reference

template <typename T>
void
sectionOfMarchingCubes(
    unsigned const& xbeg, unsigned const& ybeg, unsigned const& zbeg,
    unsigned const& xend, unsigned const& yend, unsigned const& zend,
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
    // normals and triangles. Use pointMap to not add any duplicates on
    // this thread to points or normals

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
TriangleMesh<T>
MarchingCubes(util::Image3D<T> const& image, T const& isoval, unsigned const& grainDim)
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
    std::vector<std::array<unsigned, 3> > indexTriangles;

    // If multiple threads or processes are present, points may be added
    // more than once. To fix this, duplicateTracker will be filled with pairs
    // containing the point index in points/normals and the global edge index.
    // Later, util::duplicateRemover will remove the duplicates and return
    // the final, duplicate free mesh.
    std::vector<std::pair<unsigned, unsigned> > duplicateTracker;

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

    unsigned xBeginIdx = image.xBeginIdx();
    unsigned yBeginIdx = image.yBeginIdx();
    unsigned zBeginIdx = image.zBeginIdx();

    unsigned xEndIdxExtent = image.xEndIdx();
    unsigned yEndIdxExtent = image.yEndIdx();
    unsigned zEndIdxExtent = image.zEndIdx();

    unsigned numSectX = (xEndIdxExtent - xBeginIdx + grainDim - 1) / grainDim;
    unsigned numSectY = (yEndIdxExtent - yBeginIdx + grainDim - 1) / grainDim;
    unsigned numSectZ = (zEndIdxExtent - zBeginIdx + grainDim - 1) / grainDim;

    unsigned numSections = numSectX * numSectY * numSectZ;
    unsigned numSectionsPerPage = numSectX * numSectY;

    #pragma omp parallel
    {
        // Each openMP thread manages it's own threadPoints, threadNormals,
        // threadIndexTriangles and threadPointMap.
        std::vector<std::array<T, 3> > threadPoints;
        std::vector<std::array<T, 3> > threadNormals;
        std::vector<std::array<unsigned, 3> > threadIndexTriangles;

        // threadPointMap will be a dictionary from global edge indices
        // to indices in threadPoints and threadNormals. Note that
        // threadPoints and threadNormals share the same indices
        // but are only unique among this section/thread of execution.
        std::unordered_map<unsigned, unsigned> threadPointMap;

        #pragma omp for nowait
        for(unsigned i = 0; i < numSections; ++i)
        {
            // Determine the coordinates of this section.
            unsigned xSectIdx = (i % numSectionsPerPage) % numSectX;
            unsigned ySectIdx = (i % numSectionsPerPage) / numSectX;
            unsigned zSectIdx = (i / numSectionsPerPage);

            unsigned xbeg = xSectIdx * grainDim;
            unsigned ybeg = ySectIdx * grainDim;
            unsigned zbeg = zSectIdx * grainDim;

            unsigned xend = std::min(xbeg + grainDim, xEndIdxExtent);
            unsigned yend = std::min(ybeg + grainDim, yEndIdxExtent);
            unsigned zend = std::min(zbeg + grainDim, zEndIdxExtent);

            // How does this work? TODO
            // For performance reasons, rehashing the unordered map.
            unsigned approxNumberOfEdges = 3*(xend-xbeg)*(yend-ybeg)*(zend-zbeg);
            unsigned mapSize = approxNumberOfEdges / 8 + 6;
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

            unsigned offset = points.size();

            points.insert(
                points.end(), threadPoints.begin(), threadPoints.end());
            normals.insert(
                normals.end(), threadNormals.begin(), threadNormals.end());

            // The points refered to in each tri need to refer to the points
            // in the points vector, not threadPoints.
            for(std::array<unsigned, 3>& tri: threadIndexTriangles)
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
                unsigned const& pointIndex = offset + iter->second;
                unsigned const& globalEdgeIndex = iter->first;

                duplicateTracker[pointIndex] =
                    std::make_pair(pointIndex, globalEdgeIndex);
            }
        }
    }

    // If there are no triangles whatsover, return here.
    if(indexTriangles.size() == 0)
    {
        return TriangleMesh<T>();
    }

    return util::duplicateRemover(
        duplicateTracker, points, normals, indexTriangles);
}

int main(int argc, char* argv[])
{
    char* vtkFile = argv[1];
    char* outFile = argv[2];
    float isoval = atof(argv[3]);

    // To control the granularity of the parallel execution, grainDim is passed
    // to the algorithm. grainDim is the largest number of cubes to be
    // processed in each dimension.
    unsigned grainDim = 256;
    if(argc == 5)
    {
        grainDim = atoi(argv[4]);
    }

    // Load the image file
    util::Image3D<float> image = util::loadImage<float>(vtkFile);

    // Time the output
    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();

    // MarchingCubes runs the algorithm. As inputs it takes the image
    // loaded at vtkFile and the isoval of the surface to approximate. It's
    // output is a TriangleMesh which stores the mesh as a vector
    // of triangles.
    TriangleMesh<float> polygonalMesh = MarchingCubes(image, isoval, grainDim);

    // End timing
    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();

    // Print the output
    std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
              << (c_end-c_start) / (1.0 * CLOCKS_PER_SEC) << " s\n"
              << "Wall clock time passed: "
              << std::chrono::duration<double>(t_end-t_start).count()
              << " s\n";

    saveTriangleMesh(polygonalMesh, outFile);
}
