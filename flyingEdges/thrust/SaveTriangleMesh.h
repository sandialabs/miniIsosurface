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

#include "config.h"

#include "../util/ConvertBuffer.h"
#include "../util/TypeInfo.h"

namespace util {

void
saveTriangleMesh(
    const char* fileName,
    host_vector<scalar_t> const& px,
    host_vector<scalar_t> const& py,
    host_vector<scalar_t> const& pz,
    host_vector<scalar_t> const& nx,
    host_vector<scalar_t> const& ny,
    host_vector<scalar_t> const& nz,
    host_vector<int> const& t0,
    host_vector<int> const& t1,
    host_vector<int> const& t2)
{
    std::ofstream stream(fileName);

    size_t nverts = px.size();
    size_t ntriangles = t0.size();
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

    for(int idx = 0; idx != nverts; ++idx)
    {
        *bufPointer = px[idx];
        flipEndianness(*bufPointer++);
        *bufPointer = py[idx];
        flipEndianness(*bufPointer++);
        *bufPointer = pz[idx];
        flipEndianness(*bufPointer++);
    }
    stream.write(&wbuff[0], wbuff.size());
    stream << std::endl;

    // Writing triangle indices
    bufsize = ntriangles * 4 * sizeof(size_t);
    wbuff.resize(bufsize);
    size_t *ind = reinterpret_cast<size_t*>(&wbuff[0]);

    for(int idx = 0; idx != ntriangles; ++idx)
    {
        *ind = 3;
        flipEndianness(*ind++);

        *ind = t0[idx];
        flipEndianness(*ind++);
        *ind = t1[idx];
        flipEndianness(*ind++);
        *ind = t1[idx];
        flipEndianness(*ind++);
    }

    stream << "POLYGONS " << ntriangles << " " << ntriangles * 4 << std::endl;
    stream.write(&wbuff[0], wbuff.size());
    stream << std::endl;

    // Writing normals
    bufsize = nverts * spacialDimensions * sizeof(scalar_t);
    wbuff.resize(bufsize);
    bufPointer = reinterpret_cast<scalar_t*>(&wbuff[0]);

    for(int idx = 0; idx != nverts; ++idx)
    {
        *bufPointer = nx[idx];
        flipEndianness(*bufPointer++);
        *bufPointer = ny[idx];
        flipEndianness(*bufPointer++);
        *bufPointer = nz[idx];
        flipEndianness(*bufPointer++);
    }

    stream << "POINT_DATA " << nverts << std::endl;
    stream << "NORMALS Normals " << ti.name() << std::endl;
    stream.write(&wbuff[0], wbuff.size());
    stream << std::endl;

    stream.close();
}

} // util namespace

#endif
