/*
 * loadImage.h
 *
 *  Created on: Jan 16, 2017
 *      Author: dbourge, sjmunn
 */

#ifndef IOIOIO_H_
#define IOIOIO_H_

#include <array>
#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Image3D.h"
#include "TypeInfo.h"
#include "ConvertBuffer.h"

namespace util {

template <typename T>
void
loadHeader(
    std::ifstream& stream,
    std::array<unsigned, 3>& dim,
    std::array<T, 3>& spacing,
    std::array<T, 3>& origin,
    unsigned& npoints,
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
        std::stringstream lineStream = std::stringstream(line);
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
            lineStream >> origin[0] >> origin[1] >> origin[2];
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
        std::stringstream lineStream = std::stringstream(line);
        lineStream >> tag >> npoints;
        if (tag != "POINT_DATA" || lineStream.bad())
        {
            throw bad_format("Expecting POINT_DATA <npoints>");
        }
        //if(LOG::ReportingLevel() == logDEBUG_Step)
        //{
        //    if(npoints>1000)
        //    {
        //        throw file_too_large(npoints);
        //    }
        //}
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

template <typename T>
Image3D<T>
loadImage(const char* file)
{
    std::ifstream stream(file);
    if (!stream)
        throw file_not_found(file);

    std::array<unsigned, 3> dim;
    std::array<T, 3> spacing;
    std::array<T, 3> origin;
    unsigned npoints;

    TypeInfo ti;

    // These variables are all taken by reference
    loadHeader(stream, dim, spacing, origin, npoints, ti);

    std::size_t bufsize = npoints * ti.size();
    std::vector<char> rbuf(bufsize);
    stream.read(rbuf.data(), bufsize);

    std::vector<T> data(dim[0] * dim[1] * dim[2]);

    // Converts elements in the charachter vector to elements of type
    // T and puts them into data.
    convertBufferWithTypeInfo(rbuf.data(), ti, npoints, data.data());

    stream.close();

    return Image3D<T>(data, dim, spacing, origin);
}

} // util namespace
#endif
