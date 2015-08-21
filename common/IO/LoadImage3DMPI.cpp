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
	blockExtentSet=false;
	imageDataIdx=0;
}

template<typename T>
LoadImage3DMPI<T>::LoadImage3DMPI(LoadImage3DMPI& loaderObject) {
	const unsigned * dims=loaderObject.getVolumeDimensions();

	if (dims[0]==0) throw zero_dimensions("Zero dimension error");

	xdimFile=dims[0];
	ydimFile=dims[1];
	zdimFile=dims[2];
	fileNpoints=loaderObject.getnVolumePoints();
	nPointsInBlock=0;
	typeInfo=new TypeInfo(loaderObject.typeInfo->getId());

	// Copy the line reader locally for this object
	vtkFile=loaderObject.whichFile();
	stream.open(vtkFile);

	// Skip all the header lines
	for(int iLine=0; iLine<N_HEADER_LINES;++iLine) {
		stream.ignore ( 256, '\n' );
	}
	reader = new LineStream(stream);
	blockExtentSet=false;
	imageDataIdx=0;
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
	if (blkExt[0] > xdimFile || blkExt[1] > xdimFile) throw impossible_extent("x-axis extent doesn't make sense");
	if (blkExt[2] > ydimFile || blkExt[3] > ydimFile) throw impossible_extent("y-axis extent doesn't make sense");
	if (blkExt[4] > zdimFile || blkExt[5] > zdimFile) throw impossible_extent("z-axis extent doesn't make sense");

	for (int iExt=0;iExt<6;++iExt) {
		blockExtent[iExt]=blkExt[iExt];
	}
	// add one because the march extent stops 1 unit before the end
	nPointsInBlock=blockExtent[1]-blockExtent[0]+2;
	nPointsInBlock*=blockExtent[3]-blockExtent[2]+2;
	nPointsInBlock*=blockExtent[5]-blockExtent[4]+2;
	blockExtentSet=true;
}

template<typename T>
void LoadImage3DMPI<T>::readEntireVolumeData(Image3D<T>& image) {
	unsigned blkExt[6];
	blkExt[0]=blkExt[2]=blkExt[4]=0;
	blkExt[1]=xdimFile-1;
	blkExt[3]=ydimFile-1;
	blkExt[5]=zdimFile-1;
	this->setBlockExtent(blkExt);
	this->readBlockData(image);
}

template<typename T>
void LoadImage3DMPI<T>::readBlockData(Image3D<T>& image) {
	if (!blockExtentSet) throw block_extent_not_set("Set block extent first");

	// Get to initial position
	unsigned npointsIgnore = blockExtent[0]+ (blockExtent[2] * xdimFile) + (blockExtent[4] * xdimFile*ydimFile);
	this->streamIgnore(npointsIgnore);

	unsigned nXpoints = blockExtent[1]-blockExtent[0]+2; // + 2 for first and last pts
	unsigned nYpoints = blockExtent[3]-blockExtent[2]+2;
	unsigned nZpoints = blockExtent[5]-blockExtent[4]+2;

	unsigned nXpointsIgnore = xdimFile - nXpoints;
	unsigned nYpointsIgnore = xdimFile*(ydimFile-nYpoints);

	size_t readXlineSize = nXpoints * typeInfo->size();
	size_t totalReadSize = nPointsInBlock * typeInfo->size();
	std::vector<char> rbufRead(totalReadSize);

	unsigned iYline;
	unsigned iZline;

	std::cout << "npointsIgnore " << npointsIgnore << std::endl;

	for (iZline=0;iZline < nZpoints;++iZline) {
		for (iYline=0;iYline < nYpoints;++iYline) {
			// Read a line along the x-dimension
			stream.read(&rbufRead[imageDataIdx], readXlineSize);

			// Ignore the rest of the line
			this->streamIgnore(nXpointsIgnore);
			// Update rbufRead location
			imageDataIdx+=readXlineSize;
		}
		// Ignore the x-axis lines outside the extent
		this->streamIgnore(nYpointsIgnore);
		std::cout << "Total y-ignore: " << nYpointsIgnore << std::endl;
	}

	image.setDimension(nXpoints, nYpoints, nZpoints);
	image.setSpacing(spacing[0], spacing[1], spacing[2]);
	image.setOrigin(origin[0], origin[1], origin[2]);
	image.allocate();

	convertBufferWithTypeInfo(&rbufRead[0], *typeInfo, nPointsInBlock, image.getData());
}

template<typename T>
void LoadImage3DMPI<T>::streamIgnore(unsigned nPointsIgnore) {
	size_t ignoreSize = nPointsIgnore * typeInfo->size();
	std::vector<char> rbufIgnore(ignoreSize);
	stream.read(&rbufIgnore[0], ignoreSize);
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
const unsigned& LoadImage3DMPI<T>::getnVolumePoints(void) const {
	return fileNpoints;
}

// Must instantiate class for separate compilation
template class LoadImage3DMPI<float_t> ;
