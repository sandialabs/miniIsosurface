/*
 * LoadImage3D.h
 *
 *  Created on: Jul 6, 2015
 *      Author: sjmunn
 */

#ifndef LOADIMAGE3D_H_
#define LOADIMAGE3D_H_

#include"../includes.h"

// I/O headers
#include"../Data_Obj/Image3D.h"
#include"LineStream.h"
#include"ConvertBuffer.h"

// Reporting Headers
#include"../Reporting/IO_errors.h"
#include"../Reporting/Log.h"

template<typename T>
void loadImage3D(const char *vtkFileName, Image3D<T> *image) {
	std::ifstream stream;
	stream.open(vtkFileName);
	if (!stream) {
		throw file_not_found(vtkFileName);
	}

	LineStream reader(stream);

	// read and discard header
	reader.readline();

	reader.readline();
	std::string name;
	reader.stream() >> name;

	reader.readline();
	std::string format;
	reader.stream() >> format;
	if (format != "BINARY") {
		throw bad_format("Only 'BINARY' format supported");
	}

	std::string tag;

	reader.readline();
	std::string dataset;
	reader.stream() >> tag >> dataset;
	if (tag != "DATASET" || dataset != "STRUCTURED_POINTS") {
		throw bad_format("Expecting STRUCTURED_POINTS dataset");
	}

	int count = 3;
	unsigned xdim, ydim, zdim;
	T spacing[3];
	T origin[3];
	while (count) {
		reader.readline();
		reader.stream() >> tag;

		if (tag == "DIMENSIONS") {
			reader.stream() >> xdim >> ydim >> zdim;
			if (reader.stream().bad()) {
				throw bad_format("Expecting DIMENSIONS [3]");
			}
			--count;
		} else if (tag == "SPACING") {
			reader.stream() >> spacing[0] >> spacing[1] >> spacing[2];
			if (reader.stream().bad()) {
				throw bad_format("Expecting SPACING [3]");
			}
			--count;
		} else if (tag == "ORIGIN") {
			reader.stream() >> origin[0] >> origin[1] >> origin[2];
			if (reader.stream().bad()) {
				throw bad_format("Expecting ORIGIN [3]");
			}
			--count;
		} else {
			throw bad_format("Expecting DIMENSIONS, SPACING and ORIGIN");
		}
	}

	reader.readline();
	unsigned npoints = 0;
	reader.stream() >> tag >> npoints;
	if (tag != "POINT_DATA" || reader.stream().bad()) {
		throw bad_format("Expecting POINT_DATA <npoints>");
	}

	if(LOG::ReportingLevel() == logDEBUG_Step) {
		if(npoints>1000) {
			throw file_too_large(npoints);
		}
	}

	reader.readline();
	std::string scalName, typeName;
	reader.stream() >> tag >> scalName >> typeName;
	if (tag != "SCALARS" || reader.stream().bad()) {
		throw bad_format("Expecting SCALARS <name> <type>");
	}
	TypeInfo ti = createTemplateTypeInfo(typeName.c_str());
	if (ti.getId() == TypeInfo::ID_UNKNOWN) {
		throw bad_format("Unsupported data type");
	}

	reader.readline();
	reader.stream() >> tag >> name;
	if (tag != "LOOKUP_TABLE" || reader.stream().bad()) {
		throw bad_format("Expecting LOOKUP_TABLE name");
	}
	if (name != "default") {
		unsigned size;
		reader.stream() >> size;
	}

	size_t bufsize = npoints * ti.size();
	std::vector<char> rbuf(bufsize);
	stream.read(&rbuf[0], bufsize);

	image->setDimension(xdim, ydim, zdim);
	image->setSpacing(spacing[0], spacing[1], spacing[2]);
	image->setOrigin(origin[0], origin[1], origin[2]);
	image->allocate();

	convertBufferWithTypeInfo(&rbuf[0], ti, npoints, image->getData());

	stream.close();
}

#endif /* LOADIMAGE3D_H_ */
