#include <array>
#include <vector>
#include <unordered_map>

#include <ctime>
#include <chrono>
#include <iomanip>

#include "../util/MarchingCubesTables.h"
#include "LoadImage.h"

#include "Image3D.h"
#include "TriangleMesh.h"

template<typename T> // TODO put somewhere else?
int
findCaseId(std::array<T, 8> const& cubeVertexVals, T const& isoval)
{
    int caseId=0;
    for(int i = 0; i < 8; ++i)
    {
        // Note that util::caseMask[i] = {1, 2, 4, 8, 16, 32, 64, 128}
        // Example:
        //   Suppose caseId = 3, i = 4 and cubeVertexVals[4] >= isoval.
        //   Then caseId will be set to caseId | 16. In terms of bits,
        //   caseId = 11000000 | 00001000 = 11001000 = 19
        if(cubeVertexVals[i] >= isoval)
        {
            caseId |= util::caseMask[i];
        }
    }
    return caseId;
}

template <typename T, std::size_t N> // TODO put somewhere else?
std::array<T, N>
interpolate(std::array<T, N> const& a, std::array<T, N> const& b, T weight)
{
    std::array<T, N> ret;
    for(int i = 0; i != N; ++i)
    {
        ret[i] = a[i] + (weight * (b[i] - a[i]));
    }
    return ret;
}

template <typename T>
TriangleMesh<T>
MarchingCubes(Image3D<T> const& image, T const& isoval)
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

    // A general cube has 16 edges, each labeled with a cube edge index from
    // 0 to 15. Among the whole image, each individual edge has a global
    // edge index, which can be queried from the Image3D data structure.

    // pointMap will be a dictionary from global edge indices to indices
    // in points and normals. Note that points and normals share the
    // indices.
    using PointMap = std::unordered_map<unsigned, unsigned>;
    PointMap pointMap;
    unsigned ptIdx = 0;

    unsigned xBeginIdx = image.xBeginIdx();
    unsigned yBeginIdx = image.yBeginIdx();
    unsigned zBeginIdx = image.zBeginIdx();

    // The algorithm looks at a cube at a time, where the (i, i, i)th
    // cube will contain the (i+1, i+1, i+1) vertex. So the last index
    // needed is 1 minus the end index.
    unsigned xEndIdxExtent = image.xEndIdx() - 1;
    unsigned yEndIdxExtent = image.yEndIdx() - 1;
    unsigned zEndIdxExtent = image.zEndIdx() - 1;

    // For each cube, determine whether or not the isosurface intersects
    // the given cube. If so, first find the cube configuration from a lookup
    // table. Then add the triangles of that cube configuration to points,
    // normals and triangles. Use pointMap to not add any duplicate points or
    // normals.
    for (unsigned zidx = zBeginIdx; zidx != zEndIdxExtent; ++zidx)
    {
        for (unsigned yidx = yBeginIdx; yidx != yEndIdxExtent; ++yidx)
        {
            // A buffer is used to improve cache efficency when retrieving
            // vertex values of each cube.
            auto buffer = image.createBuffer(xBeginIdx, yidx, zidx);

            for (unsigned xidx = xBeginIdx; xidx != xEndIdxExtent; ++xidx)
            {
                // For each x, y, z index, get the corresponding 8 scalar values
                // of a cube with one corner being at the x, y, z index.
                std::array<T, 8> cubeVertexVals = buffer.getCubeVertexValues(xidx);

                int cellCaseId = findCaseId(cubeVertexVals, isoval);

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
                                interpolate(posCube[v1], posCube[v2], w);
                            std::array<T, 3> newNorm =
                                interpolate(gradCube[v1], gradCube[v2], w);

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
    std::cout << indexTriangles.size() << std::endl; // TODO

    // points, normals and indexTriangles contain all the information
    // needed with respect to a new polygonal mesh, stored in the
    // TriangleMesh data structure.
    return TriangleMesh<T>(points, normals, indexTriangles);
}


int main(int argc, char* argv[])
{
    char* vtkFile = argv[1]; // TODO not implementing separate data yet
    char* outFile = argv[2];
    float isoval = atof(argv[3]); //TODO show usage error

    // Load the image file
    Image3D<float> image = loadImage<float>(vtkFile);

    // Time the output
    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();

    // MarchingCubes runs the algorithm. As inputs it takes the image
    // loaded at vtkFile and the isoval of the surface to approximate. It's
    // output is a TriangleMesh which stores the mesh as a vector
    // of triangles.
    TriangleMesh<float> polygonalMesh = MarchingCubes(image, isoval);

    // End timing
    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();

    // Print the output
    std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
              << (c_end-c_start) / (1.0 * CLOCKS_PER_SEC) << " s\n"
              << "Wall clock time passed: "
              << std::chrono::duration<double>(t_end-t_start).count()
              << " s\n";

    //saveTriangleMesh(mesh, outFile); TODO
}
