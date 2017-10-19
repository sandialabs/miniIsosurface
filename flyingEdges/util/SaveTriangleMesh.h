/*
 * SaveTriangleMesh.h
 *
 *  Created on: Jan 18, 2017
 *      Author: dbourge, sjmunn
 *
 * miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
 * See LICENSE.txt for details.
 *
 * Copyright (c) 2017
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
 * the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
 * certain rights in this software.
 */

#ifndef UTIL_SAVETRIANGLEMESH_H_
#define UTIL_SAVETRIANGLEMESH_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "FlyingEdges_Config.h"

#include "ConvertBuffer.h"
#include "TypeInfo.h"

#include "TriangleMesh.h"

namespace util {

void
saveTriangleMesh(TriangleMesh const& mesh, const char* fileName)
{
    using PointIterator =    typename TriangleMesh::PointIterator;
    using TriangleIterator = typename TriangleMesh::TriangleIterator;
    using NormalIterator =   typename TriangleMesh::NormalIterator;

    std::ofstream stream(fileName);

    size_t nverts = mesh.numberOfVertices();
    size_t ntriangles = mesh.numberOfTriangles();
    size_t spacialDimensions = 3;

    TypeInfo ti = createTemplateTypeInfo<scalar_t>();

    stream << "# vtk DataFile Version 3.0" << std::endl;
    stream << "Isosurface Mesh" << std::endl;
    stream << "BINARY" << std::endl;
    stream << "DATASET POLYDATA" << std::endl;
    stream << "POINTS " << nverts << " " << ti.name() << std::endl;

    std::vector<char> wbuff;
    std::size_t bufsize = nverts * spacialDimensions * sizeof(scalar_t);
    wbuff.resize(bufsize);

    // Writing points data
    scalar_t *bufPointer = reinterpret_cast<scalar_t*>(&wbuff[0]);

    PointIterator ptBegIter = mesh.pointsBegin();
    PointIterator ptEndIter = mesh.pointsEnd();
    for(PointIterator iter = ptBegIter; iter != ptEndIter; ++iter)
    {
        for (int i = 0; i < 3; ++i)
        {
            *bufPointer = (*iter)[i];
            flipEndianness(*bufPointer++);
        }
    }
    stream.write(&wbuff[0], wbuff.size());
    stream << std::endl;

    // Writing triangle indices
    bufsize = ntriangles * 4 * sizeof(size_t);
    wbuff.resize(bufsize);
    size_t *ind = reinterpret_cast<size_t*>(&wbuff[0]);

    TriangleIterator triBegIter = mesh.trianglesBegin();
    TriangleIterator triEndIter = mesh.trianglesEnd();
    for(TriangleIterator iter = triBegIter; iter != triEndIter; ++iter)
    {
        *ind = 3;
        flipEndianness(*ind++);
        for (int i = 0; i < 3; ++i)
        {
            *ind = (*iter)[i];
            flipEndianness(*ind++);
        }
    }

    stream << "POLYGONS " << ntriangles << " " << ntriangles * 4 << std::endl;
    stream.write(&wbuff[0], wbuff.size());
    stream << std::endl;

    // Writing normals
    bufsize = nverts * spacialDimensions * sizeof(scalar_t);
    wbuff.resize(bufsize);
    bufPointer = reinterpret_cast<scalar_t*>(&wbuff[0]);

    NormalIterator norBegIter = mesh.normalsBegin();
    NormalIterator norEndIter = mesh.normalsEnd();
    for(NormalIterator iter = norBegIter; iter != norEndIter; ++iter)
    {
        for (int i = 0; i < 3; ++i)
        {
            *bufPointer = (*iter)[i];
            flipEndianness(*bufPointer++);
        }
    }

    stream << "POINT_DATA " << nverts << std::endl;
    stream << "NORMALS Normals " << ti.name() << std::endl;
    stream.write(&wbuff[0], wbuff.size());
    stream << std::endl;

    stream.close();
}

} // util namespace

#endif
