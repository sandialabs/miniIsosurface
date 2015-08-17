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
void SerialAlgo<T>::visit(SerialData<T> &data){
	std::cout << "Hello World" << std::endl;

	const unsigned *dims = data.imageIn.getDimension();
	const T *origin = data.imageIn.getOrigin();
	const T *spacing = data.imageIn.getSpacing();

	CLOG(logYAML) << "Marching cubes algorithm: SERIAL";
	data.doc.add("Marching cubes algorithm", "SERIAL");

	T range[6];
	for (int i = 0; i < 3; ++i) {
		range[i * 2] = origin[i];
		range[(i * 2) + 1] = origin[i]
				+ (static_cast<T>(dims[i] - 1) * spacing[i]);
	}

	Range3D cellRange(0, dims[2] - 1, grainDim, 0, dims[1] - 1, grainDim, 0,
			dims[0] - 1, grainDim);
	cellRange.extent(data.ext);

	data.initEdgeIndices();

	unsigned mapSize = data.edgeIndices->nAllEdges / 8; // Very approximate hack..
	data.pointMap.rehash(mapSize);

	MarchAlgorithm<T>::extractIsosurfaceFromBlock(data.imageIn, data.ext, data.isoval, data.pointMap, data.edgeIndices, &data.mesh);
}


// Must instantiate class for separate compilation
template class SerialAlgo<float_t> ;
