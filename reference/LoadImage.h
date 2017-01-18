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
#include "../util/TypeInfo.h"
#include "../util/ConvertBuffer.h"

template <typename T>
Image3D<T>
loadImage(const char* file)
{
    std::array<unsigned, 3> dim;
    std::array<unsigned, 3> spacing;
    std::array<unsigned, 3> origin;

    std::ifstream stream(file);
    if (!stream)
        throw util::file_not_found(file);

    std::string line;
    std::string tag;

    // Discarding these lines
    std::getline(stream, line);
    std::getline(stream, line);

    std::getline(stream, line);
    std::string format = line;
    if (format != "BINARY")
    {
        throw util::bad_format("Only 'BINARY' format supported");
    }

    {
        std::getline(stream, line);
        std::stringstream lineStream(line);
        std::string dataset;
        lineStream >> tag >> dataset;
        if (tag != "DATASET" || dataset != "STRUCTURED_POINTS" || lineStream.bad())
        {
            throw util::bad_format("Expecting STRUCTURED_POINTS dataset");
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
                throw util::bad_format("Expecting DIMENSIONS [3]");
            }
        }
        else if (tag == "SPACING")
        {
            lineStream >> spacing[0] >> spacing[1] >> spacing[2];
            if (lineStream.bad())
            {
                throw util::bad_format("Expecting SPACING [3]");
            }
        }
        else if (tag == "ORIGIN")
        {
            lineStream >> origin[0] >> origin[1] >> origin[2];
            if (lineStream.bad())
            {
                throw util::bad_format("Expecting ORIGIN [3]");
            }
        }
        else
        {
            throw util::bad_format("Expecting DIMENSIONS, SPACING and ORIGIN");
        }
    }

    unsigned npoints = 0;
    {
        std::getline(stream, line);
        std::stringstream lineStream = std::stringstream(line);
        lineStream >> tag >> npoints;
        if (tag != "POINT_DATA" || lineStream.bad())
        {
            throw util::bad_format("Expecting POINT_DATA <npoints>");
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
            throw util::bad_format("Expecting SCALARS <name> <type>");
        }
    }

    util::TypeInfo ti(typeName.c_str());
    if (ti.getId() == util::TypeInfo::ID_UNKNOWN)
    {
        throw util::bad_format("Unsupported data type");
    }

    {
        std::getline(stream, line);
        std::stringstream lineStream(line);
        std::string name;
        lineStream >> tag >> name;
        if (tag != "LOOKUP_TABLE" || lineStream.bad())
        {
            throw util::bad_format("Expecting LOOKUP_TABLE name");
        }
    }

    std::size_t bufsize = npoints * ti.size();
    std::vector<char> rbuf(bufsize);
    stream.read(rbuf.data(), bufsize);

    std::vector<T> data(dim[0] * dim[1] * dim[2]);

    // Converts elements in the charachter vector to elements of type
    // T and puts them into data.
    util::convertBufferWithTypeInfo(rbuf.data(), ti, npoints, data.data());

    stream.close();

    return Image3D<T>(data, dim, spacing, origin);
}

#endif
