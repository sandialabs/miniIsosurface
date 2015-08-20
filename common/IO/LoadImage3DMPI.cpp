/*
 * LoadImage3DMPI.cpp
 *
 *  Created on: Aug 20, 2015
 *      Author: sjmunn
 */

#include "LoadImage3DMPI.h"

template<typename T>
LoadImage3DMPI<T>::LoadImage3DMPI() {
	xdimFile=0;
	ydimFile=0;
	zdimFile=0;
	fileNpoints=0;
	reader=0;
	vtkFile=0;
	nPointsInBlock=0;
	typeInfo=0;
}

template<typename T>
LoadImage3DMPI<T>::LoadImage3DMPI(LoadImage3DMPI& loaderObject) {
	const unsigned * dims=loaderObject.getVolumeDimensions();
	xdimFile=dims[0];
	ydimFile=dims[1];
	zdimFile=dims[2];
	fileNpoints=loaderObject.getnVolumePoints();
	nPointsInBlock=0;
	BOOST_TEST_CHECKPOINT("Copied header info");
	typeInfo=new TypeInfo(loaderObject.typeInfo->getId());

	BOOST_TEST_CHECKPOINT("Copied type info");
	// Copy the line reader locally for this object
	vtkFile=loaderObject.whichFile();
	stream.open(vtkFile);

	BOOST_TEST_CHECKPOINT("Opened file, skipping hear data");
	// Skip all the header lines
	for(int iLine=0; iLine<N_HEADER_LINES;++iLine) {
		stream.ignore ( std::numeric_limits<std::streamsize>::max(), '\n' );
	}
	BOOST_TEST_CHECKPOINT("Initializing line stream at binary data");
	reader = new LineStream(stream);
}

template<typename T>
LoadImage3DMPI<T>::~LoadImage3DMPI() {
	stream.close();
	delete reader;
	delete typeInfo;
}

template<typename T>
void LoadImage3DMPI<T>::loadHeader(const char *vtkFileName) {
	vtkFile=vtkFileName;
	stream.open(vtkFileName);
	if (!stream) {
		throw file_not_found(vtkFileName);
	}

	reader = new LineStream(stream);

	// read and discard header
	reader->readline();

	reader->readline();
	std::string name;
	reader->stream() >> name;

	reader->readline();
	std::string format;
	reader->stream() >> format;
	if (format != "BINARY") {
		throw bad_format("Only 'BINARY' format supported");
	}

	std::string tag;

	reader->readline();
	std::string dataset;
	reader->stream() >> tag >> dataset;
	if (tag != "DATASET" || dataset != "STRUCTURED_POINTS") {
		throw bad_format("Expecting STRUCTURED_POINTS dataset");
	}

	int count = 3;
	while (count) {
		reader->readline();
		reader->stream() >> tag;

		if (tag == "DIMENSIONS") {
			reader->stream() >> xdimFile >> ydimFile >> zdimFile;
			if (reader->stream().bad()) {
				throw bad_format("Expecting DIMENSIONS [3]");
			}
			--count;
		} else if (tag == "SPACING") {
			reader->stream() >> spacing[0] >> spacing[1] >> spacing[2];
			if (reader->stream().bad()) {
				throw bad_format("Expecting SPACING [3]");
			}
			--count;
		} else if (tag == "ORIGIN") {
			reader->stream() >> origin[0] >> origin[1] >> origin[2];
			if (reader->stream().bad()) {
				throw bad_format("Expecting ORIGIN [3]");
			}
			--count;
		} else {
			throw bad_format("Expecting DIMENSIONS, SPACING and ORIGIN");
		}
	}

	reader->readline();
	unsigned npoints = 0;
	reader->stream() >> tag >> npoints;
	if (tag != "POINT_DATA" || reader->stream().bad()) {
		throw bad_format("Expecting POINT_DATA <npoints>");
	}

	if(LOG::ReportingLevel() == logDEBUG_Step) {
		if(npoints>1000) {
			throw file_too_large(npoints);
		}
	}
	fileNpoints=npoints;

	reader->readline();
	std::string scalName, typeName;
	reader->stream() >> tag >> scalName >> typeName;
	if (tag != "SCALARS" || reader->stream().bad()) {
		throw bad_format("Expecting SCALARS <name> <type>");
	}

	TypeInfo newTypeInfo = createTypeInfo(typeName.c_str());
	typeInfo = new TypeInfo(newTypeInfo.getId()); // newTypeInfo has limited scope
	if (typeInfo->getId() == TypeInfo::ID_UNKNOWN) {
		throw bad_format("Unsupported data type");
	}

	reader->readline();
	reader->stream() >> tag >> name;
	if (tag != "LOOKUP_TABLE" || reader->stream().bad()) {
		throw bad_format("Expecting LOOKUP_TABLE name");
	}
	if (name != "default") {
		unsigned size;
		reader->stream() >> size;
	}
}

template<typename T>
void LoadImage3DMPI<T>::setBlockExtent(const unsigned * blkExt) {
	for (int iExt=0;iExt<6;++iExt) {
		blockExtent[iExt]=blkExt[iExt];
	}
	// add one because the march extent stops 1 unit before the end
	nPointsInBlock=blockExtent[1]-blockExtent[0]+1;
	nPointsInBlock+=blockExtent[3]-blockExtent[2]+1;
	nPointsInBlock+=blockExtent[5]-blockExtent[4]+1;
}

template<typename T>
void LoadImage3DMPI<T>::readBlockData(Image3D<T>& image) {
//	size_t bufsize = nPointsInBlock * typeInfo.size();
//	std::vector<char> rbuf(bufsize);
//	stream.read(&rbuf[0], bufsize);
	// Get to initial position
	unsigned npointsIgnore = blockExtent[0]+ (blockExtent[2] * xdimFile) + (blockExtent[4] * xdimFile*ydimFile);
	size_t ignoreSize = npointsIgnore * typeInfo->size();
	std::vector<char> rbufIgnore(ignoreSize);
	BOOST_TEST_CHECKPOINT("Skip data file to origin, npointsIgnore " << npointsIgnore );
	stream.read(&rbufIgnore[0], ignoreSize);

	// DEBUG just read one point
	unsigned npointsRead = 1;
	size_t readSize = npointsRead * typeInfo->size();
	std::vector<char> rbufRead(readSize);
	BOOST_TEST_CHECKPOINT("Get actual data, nPoints " << npointsRead);
	stream.read(&rbufRead[0], readSize);

	image.setDimension(1, 1, 1);
	image.setSpacing(spacing[0], spacing[1], spacing[2]);
	image.setOrigin(origin[0], origin[1], origin[2]);
	image.allocate();

	BOOST_TEST_CHECKPOINT("Converting buffer");
	convertBufferWithTypeInfo(&rbufRead[0], *typeInfo, npointsRead, image.getData());
}

template<typename T>
const unsigned* LoadImage3DMPI<T>::getVolumeDimensions(void) const {
	unsigned *dims = new unsigned[3];
	dims[0]=xdimFile;
	dims[1]=ydimFile;
	dims[2]=zdimFile;
	return &dims[0];
}

template<typename T>
const unsigned LoadImage3DMPI<T>::getnVolumePoints(void) const {
	return fileNpoints;
}


