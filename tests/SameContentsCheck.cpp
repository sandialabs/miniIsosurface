#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>
#include <array>

#include <algorithm>
#include <unordered_map>

#include "../reference/TriangleMesh.h"
#include "../util/ConvertBuffer.h"

// Quick test used to make sure that output files are exactly
// equal. Use this to make sure that the new serial version
// matches with the old version. Namely that saveTriangleMesh
// is correct.
bool checkSameFiles(char* fileA, char* fileB)
{
    std::ifstream streamA(fileA);
    std::ifstream streamB(fileB);

    std::string lineA("asd");
    std::string lineB("zxc");

    if(!streamA) return false;
    if(!streamB) return false;

    while(!streamA.eof() && !streamB.eof())
    {
        std::getline(streamA, lineA);
        std::getline(streamB, lineB);

        if (lineA != lineB)
            return false;
    }

    return streamA.eof() && streamB.eof();
}

TriangleMesh<float> LoadFloatMesh(char* file)
{
    std::vector<std::array<float, 3> > points;
    std::vector<std::array<unsigned, 3> > indexTriangles;
    std::vector<std::array<float, 3> > normals;

    std::ifstream stream(file);

    std::string line;
    std::string discard;

    std::getline(stream, line);
    std::getline(stream, line);
    std::getline(stream, line);
    std::getline(stream, line);

    unsigned nverts;
    {
        std::getline(stream, line);
        std::stringstream lineStream(line);

        lineStream >> discard >> nverts >> discard;
    }
    points.resize(nverts);
    normals.resize(nverts);

    std::vector<char> wbuff;
    unsigned spatialDimensions = 3;
    unsigned floatSize = 4;
    std::size_t bufsize = nverts * spatialDimensions * floatSize;
    wbuff.resize(bufsize);
    {
        //std::getline(stream, line);
        //std::stringstream lineStream(line);
        //lineStream.read(&wbuff[0], wbuff.size());
        stream.read(&wbuff[0], wbuff.size());

        float* bufPointer = reinterpret_cast<float*>(&wbuff[0]);
        for(std::array<float, 3>& point: points)
        {
            util::flipEndianness(bufPointer[0]);
            util::flipEndianness(bufPointer[1]);
            util::flipEndianness(bufPointer[2]);

            point[0] = bufPointer[0];
            point[1] = bufPointer[1];
            point[2] = bufPointer[2];

            bufPointer += 3;
        }

        std::getline(stream, line);
    }

    unsigned ntriangles;
    {
        std::getline(stream, line);
        std::stringstream lineStream(line);

        lineStream >> discard;
        lineStream >> ntriangles;
        lineStream >> discard;

        lineStream >> discard >> ntriangles >> discard;
    }
    indexTriangles.resize(ntriangles);

    bufsize = ntriangles * 4 * sizeof(unsigned);
    wbuff.resize(bufsize);
    {
        stream.read(&wbuff[0], wbuff.size());

        unsigned* bufPointer = reinterpret_cast<unsigned*>(&wbuff[0]);
        for(std::array<unsigned, 3>& tri: indexTriangles)
        {
            util::flipEndianness(bufPointer[1]);
            util::flipEndianness(bufPointer[2]);
            util::flipEndianness(bufPointer[3]);

            tri[0] = bufPointer[1];
            tri[1] = bufPointer[2];
            tri[2] = bufPointer[3];

            bufPointer += 4;
        }

        std::getline(stream, line);
    }

    std::getline(stream, line);
    std::getline(stream, line);

    bufsize = nverts * spatialDimensions * floatSize;
    wbuff.resize(bufsize);
    {
        stream.read(&wbuff[0], wbuff.size());

        float* bufPointer = reinterpret_cast<float*>(&wbuff[0]);
        for(std::array<float, 3>& normal: normals)
        {
            util::flipEndianness(bufPointer[0]);
            util::flipEndianness(bufPointer[1]);
            util::flipEndianness(bufPointer[2]);

            normal[0] = bufPointer[0];
            normal[1] = bufPointer[1];
            normal[2] = bufPointer[2];

            bufPointer += 3;
        }
    }

    return TriangleMesh<float>(points, normals, indexTriangles);
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
    std::unordered_map<std::array<T, 3>, unsigned, arrayHash<T, 3> >;

bool sameMesh(TriangleMesh<float> const& meshA, TriangleMesh<float> const& meshB)
{
    if(meshA.numberOfTriangles() != meshB.numberOfTriangles())
    {
        std::cout << meshA.numberOfTriangles() << ", " << meshB.numberOfTriangles() << std::endl;
        std::cout << "not the same number of triangles" << std::endl;
        return false;
    }

    unsigned numTris = meshA.numberOfTriangles();


    // They may not have the same number of vertices!
    /*
    if(meshA.numberOfVertices() != meshB.numberOfVertices())
    {
        std::cout << "not the same number of vertices" << std::endl;
        return false;
    }

    unsigned d1 = std::distance(meshA.normalsBegin(), meshA.normalsEnd());
    unsigned d2 = std::distance(meshB.normalsBegin(), meshB.normalsEnd());
    if(d1 != d2)
    {
        std::cout << "not the same number of normals" << std::endl;
        return false;
    }
    */

    UnorderedMapArr<float> pointIdxA, pointIdxB;
    UnorderedMapArr<float> normalIdxA, normalIdxB;

    auto pointsA = meshA.pointsBegin();
    auto pointsB = meshB.pointsBegin();
    auto normalsA = meshA.normalsBegin();
    auto normalsB = meshB.normalsBegin();

    ///////////////////////////////////////////////////////////////////////////
    // Check that all points and normals in meshA are in meshB
    for(unsigned i = 0; i != meshB.numberOfVertices(); ++i)
    {
        pointIdxB[pointsB[i]] = i;
        normalIdxB[normalsB[i]] = i;
    }

    for(unsigned i = 0; i != meshA.numberOfVertices(); ++i)
    {
        auto itFoundPoint = pointIdxB.find(pointsA[i]);
        if(itFoundPoint == pointIdxB.end())
        {
            std::cout << "point not found" << std::endl;
            return false;
        }

        auto itFoundNormal = normalIdxB.find(normalsA[i]);
        if(itFoundNormal == normalIdxB.end())
        {
            std::cout << "normal not found" << std::endl;
            return false;
        }
    }

    // Check that all points and normals in meshB are in meshA
    for(unsigned i = 0; i != meshA.numberOfVertices(); ++i)
    {
        pointIdxA[pointsA[i]] = i;
        normalIdxA[normalsA[i]] = i;
    }

    for(unsigned i = 0; i != meshB.numberOfVertices(); ++i)
    {
        auto itFoundPoint = pointIdxA.find(pointsB[i]);
        if(itFoundPoint == pointIdxA.end())
        {
            std::cout << "point not found" << std::endl;
            return false;
        }

        auto itFoundNormal = normalIdxA.find(normalsB[i]);
        if(itFoundNormal == normalIdxA.end())
        {
            std::cout << "normal not found" << std::endl;
            return false;
        }
    }
    //
    ///////////////////////////////////////////////////////////////////////////

    auto trianglesA = meshA.trianglesBegin();
    auto trianglesB = meshB.trianglesBegin();

    std::unordered_map<std::array<float, 9>, unsigned, arrayHash<float, 9> > triIdxB;

    // The triangles don't have to have the same indexing scheme because it could be
    // the case that a point appears multiple times in points. This is a byproduct
    // of merging points and normals for parallel runs of the algorithm.

    for(unsigned i = 0; i != numTris; ++i)
    {
        float a0 = pointsB[trianglesB[i][0]][0];
        float a1 = pointsB[trianglesB[i][0]][1];
        float a2 = pointsB[trianglesB[i][0]][2];

        float a3 = pointsB[trianglesB[i][1]][0];
        float a4 = pointsB[trianglesB[i][1]][1];
        float a5 = pointsB[trianglesB[i][1]][2];

        float a6 = pointsB[trianglesB[i][2]][0];
        float a7 = pointsB[trianglesB[i][2]][1];
        float a8 = pointsB[trianglesB[i][2]][2];

        std::array<float, 9> arr = {a0, a1, a2, a3, a4, a5, a6, a7, a8};
        triIdxB[arr] = i;
    }

    for(unsigned i = 0; i != numTris; ++i)
    {
        std::array<unsigned, 3> triA = trianglesA[i];

        float a0 = pointsA[triA[0]][0];
        float a1 = pointsA[triA[0]][1];
        float a2 = pointsA[triA[0]][2];

        float a3 = pointsA[triA[1]][0];
        float a4 = pointsA[triA[1]][1];
        float a5 = pointsA[triA[1]][2];

        float a6 = pointsA[triA[2]][0];
        float a7 = pointsA[triA[2]][1];
        float a8 = pointsA[triA[2]][2];

        std::array<float, 9> arr = {a0, a1, a2, a3, a4, a5, a6, a7, a8};

        auto itFoundTri = triIdxB.find(arr);
        if(itFoundTri == triIdxB.end())
        {
            std::cout << "triangle not found" << std::endl;
            return false;
        }

    }

    return true;
}

int main(int argc, char* argv[])
{
    char* fileA = argv[1];
    char* fileB = argv[2];

    TriangleMesh<float> meshA = LoadFloatMesh(fileA);
    TriangleMesh<float> meshB = LoadFloatMesh(fileB);

    if(sameMesh(meshA, meshB))
        std::cout << "The two meshes are equivalent." << std::endl;
    else
        std::cout << "ERROR: The two meshes are not equivalent." << std::endl;
}
