/*
 * loadImage.h
 *
 *  Created on: Jan 16, 2017
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

#ifndef IOIOIO_H_
#define IOIOIO_H_

#include <array>
#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include "config.h"

#include "../util/TypeInfo.h"
#include "../util/ConvertBuffer.h"

namespace util {

void
loadHeader(
    std::ifstream& stream,
    std::array<size_t, 3>& dim,
    std::array<scalar_t, 3>& spacing,
    std::array<scalar_t, 3>& zeroPos,
    size_t& npoints,
    TypeInfo& ti)
{
    std::string line;
    std::string tag;

    // Discarding these lines
    std::getline(stream, line);
    std::getline(stream, line);

    std::getline(stream, line);
    std::string format = line;
    if (format != "BINARY")
    {
        throw bad_format("Only 'BINARY' format supported");
    }

    {
        std::getline(stream, line);
        std::stringstream lineStream(line);
        std::string dataset;
        lineStream >> tag >> dataset;
        if (tag != "DATASET" || dataset != "STRUCTURED_POINTS" || lineStream.bad())
        {
            throw bad_format("Expecting STRUCTURED_POINTS dataset");
        }
    }

    for(int count = 0; count != 3; ++count)
    {
        std::getline(stream, line);
        std::stringstream lineStream(line);
        lineStream >> tag;

        if (tag == "DIMENSIONS")
        {
            lineStream >> dim[0] >> dim[1] >> dim[2];
            if (lineStream.bad())
            {
                throw bad_format("Expecting DIMENSIONS [3]");
            }
        }
        else if (tag == "SPACING")
        {
            lineStream >> spacing[0] >> spacing[1] >> spacing[2];
            if (lineStream.bad())
            {
                throw bad_format("Expecting SPACING [3]");
            }
        }
        else if (tag == "ORIGIN")
        {
            lineStream >> zeroPos[0] >> zeroPos[1] >> zeroPos[2];
            if (lineStream.bad())
            {
                throw bad_format("Expecting ORIGIN [3]");
            }
        }
        else
        {
            throw bad_format("Expecting DIMENSIONS, SPACING and ORIGIN");
        }
    }

    {
        std::getline(stream, line);
        std::stringstream lineStream(line);
        lineStream >> tag >> npoints;
        if (tag != "POINT_DATA" || lineStream.bad())
        {
            throw bad_format("Expecting POINT_DATA <npoints>");
        }
    }

    std::string typeName;
    {
        std::getline(stream, line);
        std::stringstream lineStream(line);
        std::string scalName;
        lineStream >> tag >> scalName >> typeName;
        if (tag != "SCALARS" || lineStream.bad())
        {
            throw bad_format("Expecting SCALARS <name> <type>");
        }
    }

    ti = TypeInfo(typeName.c_str());
    if (ti.getId() == TypeInfo::ID_UNKNOWN)
    {
        throw bad_format("Unsupported data type");
    }

    {
        std::getline(stream, line);
        std::stringstream lineStream(line);
        std::string name;
        lineStream >> tag >> name;
        if (tag != "LOOKUP_TABLE" || lineStream.bad())
        {
            throw bad_format("Expecting LOOKUP_TABLE name");
        }
    }
}

void skipHeader(std::ifstream& stream)
{
    std::string line;
    for(int i = 0; i != 10; ++i)
        std::getline(stream, line);
}

void
streamIgnore(std::ifstream& stream, size_t nPoints, std::size_t pointSize)
{
    size_t increment=134217728; // Read at most 1 GB at a time

    for(size_t iIgnore = 0; iIgnore < nPoints; iIgnore += increment)
    {
        size_t diffIgnore = nPoints - iIgnore;

        if(diffIgnore > increment)
        {
            std::size_t ignoreSize = increment * pointSize;
            std::vector<char> rbufIgnore(ignoreSize);
            stream.read(&rbufIgnore[0], ignoreSize);
        }
        else
        {
            std::size_t ignoreSize = diffIgnore * pointSize;
            std::vector<char> rbufIgnore(ignoreSize);
            stream.read(&rbufIgnore[0], ignoreSize);
        }
    }
}

void loadImage_thrust(
    const char* file,
    std::vector<scalar_t>& data,
    scalar_t& spacingX,
    scalar_t& spacingY,
    scalar_t& spacingZ,
    scalar_t& zeroPosX,
    scalar_t& zeroPosY,
    scalar_t& zeroPosZ,
    int& nx,
    int& ny,
    int& nz)
{
    std::ifstream stream(file);
    if (!stream)
        throw file_not_found(file);

    std::array<size_t, 3> dim;
    std::array<scalar_t, 3> spacing;
    std::array<scalar_t, 3> zeroPos;
    size_t npoints;

    TypeInfo ti;

    // These variables are all taken by reference
    loadHeader(stream, dim, spacing, zeroPos, npoints, ti);

    std::size_t bufsize = npoints * ti.size();
    std::vector<char> rbuf(bufsize);
    stream.read(rbuf.data(), bufsize);

    data.resize(dim[0] * dim[1] * dim[2]);

    // Converts elements in the charachter vector to elements of type
    // T and puts them into data.
    convertBufferWithTypeInfo(rbuf.data(), ti, npoints, data.data());

    stream.close();

    spacingX = spacing[0];
    spacingY = spacing[1];
    spacingZ = spacing[2];
    zeroPosX = zeroPos[0];
    zeroPosY = zeroPos[1];
    zeroPosZ = zeroPos[2];
    nx = dim[0];
    ny = dim[1];
    nz = dim[2];
}

} // util namespace

#endif
