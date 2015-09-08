/*
 * SerialAlgo.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "SerialAlgo.h"

template<typename T>
SerialAlgo<T>::SerialAlgo() {
	// TODO Auto-generated constructor stub

}

template<typename T>
SerialAlgo<T>::~SerialAlgo() {
	// TODO Auto-generated destructor stub
}

template<typename T>
void SerialAlgo<T>::march(GeneralContext<T> &data) {

	const unsigned *dims = data.imageIn.getDimension();

	CLOG(logYAML) << "Marching cubes algorithm: SERIAL";
	data.doc.add("Marching cubes algorithm", "SERIAL");

	unsigned grainDim = dims[0]; // This doesn't matter

	Range3D cellRange(0, dims[2] - 1, grainDim, 0, dims[1] - 1, grainDim, 0,
			dims[0] - 1, grainDim);
	cellRange.extent(data.ext);

	setGlobalVariables(data);

	// The block is simply the full data extent
	MarchAlgorithm<T>::extractIsosurfaceFromBlock(data.imageIn, data.ext,
			data.isoval, this->globalPointMap, *(this->globalEdgeIndices), data.mesh);
}

// Must instantiate class for separate compilation
template class SerialAlgo<float_t> ;
