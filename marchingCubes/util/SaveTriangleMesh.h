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

#include "../util/ConvertBuffer.h"
#include "../util/TypeInfo.h"

#include "TriangleMesh.h"

using std::size_t;

namespace util {

template<typename T>
void
saveTriangleMesh(TriangleMesh<T> const& mesh, const char* fileName)
{
    using PointIterator =    typename TriangleMesh<T>::PointIterator;
    using TriangleIterator = typename TriangleMesh<T>::TriangleIterator;
    using NormalIterator =   typename TriangleMesh<T>::NormalIterator;

    std::ofstream stream(fileName);

    size_t nverts = mesh.numberOfVertices();
    size_t ntriangles = mesh.numberOfTriangles();
    size_t spacialDimensions = 3;

    TypeInfo ti = createTemplateTypeInfo<T>();

    stream << "# vtk DataFile Version 3.0" << std::endl;
    stream << "Isosurface Mesh" << std::endl;
    stream << "BINARY" << std::endl;
    stream << "DATASET POLYDATA" << std::endl;
    stream << "POINTS " << nverts << " " << ti.name() << std::endl;

    std::vector<char> wbuff;
    std::size_t bufsize = nverts * spacialDimensions * sizeof(T);
    wbuff.resize(bufsize);

    // Writing points data
    T *bufPointer = reinterpret_cast<T*>(&wbuff[0]);

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
    bufsize = nverts * spacialDimensions * sizeof(T);
    wbuff.resize(bufsize);
    bufPointer = reinterpret_cast<T*>(&wbuff[0]);

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
