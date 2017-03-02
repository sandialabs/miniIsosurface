/*
 * serial/main.cpp
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

#include "../util/Image3D.h"
#include "../util/TriangleMesh.h"

#include "../util/LoadImage.h"
#include "../util/SaveTriangleMesh.h"

#include "../util/util.h" // util::findCaseId, util::interpolate
#include "../util/MarchingCubesTables.h"

#include "../util/Timer.h"
#include "../mantevoCommon/YAML_Doc.hpp"

using std::size_t;

template <typename T>
util::TriangleMesh<T>
MarchingCubes(util::Image3D<T> const& image, T const& isoval)
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

    // A general cube has 16 edges, each labeled with a cube edge index from
    // 0 to 15. Among the whole image, each individual edge has a global
    // edge index, which can be queried from the Image3D data structure.

    // pointMap will be a dictionary from global edge indices to indices
    // in points and normals. Note that points and normals share the
    // indices.
    using PointMap = std::unordered_map<size_t, size_t>;
    PointMap pointMap;
    size_t ptIdx = 0;

    size_t xBeginIdx = image.xBeginIdx();
    size_t yBeginIdx = image.yBeginIdx();
    size_t zBeginIdx = image.zBeginIdx();

    // The algorithm looks at a cube ata a time, where the (i, i, i)th
    // cube will contain the (i+1, i+1, i+1) vertex. Since Image3D
    // is used to access cube information, the EndIdx returns the
    // one value past the last i such that i+1 is still within the image.
    // So if the image contains 512 points in the x dimension, the
    // valid index range is from 0 to 510. image.xEndIdx() returns one
    // past that, so image.xEndIdx() returns 511.
    size_t xEndIdx = image.xEndIdx();
    size_t yEndIdx = image.yEndIdx();
    size_t zEndIdx = image.zEndIdx();

    // For each cube, determine whether or not the isosurface intersects
    // the given cube. If so, first find the cube configuration from a lookup
    // table. Then add the triangles of that cube configuration to points,
    // normals and triangles. Use pointMap to not add any duplicate points or
    // normals.
    for (size_t zidx = zBeginIdx; zidx != zEndIdx; ++zidx)
    {
        for (size_t yidx = yBeginIdx; yidx != yEndIdx; ++yidx)
        {
            // A buffer is used to improve cache efficency when retrieving
            // vertex values of each cube.
            auto buffer = image.createBuffer(xBeginIdx, yidx, zidx);

            for (size_t xidx = xBeginIdx; xidx != xEndIdx; ++xidx)
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
    // points, normals and indexTriangles contain all the information
    // needed with respect to this new polygonal mesh, stored in TriangleMesh.
    return util::TriangleMesh<T>(points, normals, indexTriangles);
}

int main(int argc, char* argv[])
{
    float isoval;
    bool isovalSet = false;
    char* vtkFile = NULL;
    char* outFile = NULL;
    std::string yamlDirectory = "";
    std::string yamlFileName  = "";

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
    doc.add("Marching Cubes Algorithm", "serial");
    doc.add("Volume image data file path", vtkFile);
    doc.add("Polygonal mesh output file", outFile);
    doc.add("Isoval", isoval);

    // Load the image file
    util::Image3D<float> image = util::loadImage<float>(vtkFile);

    doc.add("File x-dimension", image.xdimension());
    doc.add("File y-dimension", image.ydimension());
    doc.add("File z-dimension", image.zdimension());

    // Time the output. Timer's constructor starts timing.
    util::Timer runTime;

    // MarchingCubes runs the algorithm. As inputs it takes the image
    // loaded at vtkFile and the isoval of the surface to approximate. It's
    // output is a TriangleMesh which stores the mesh as a vector
    // of triangles.
    util::TriangleMesh<float> polygonalMesh = MarchingCubes(image, isoval);

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
