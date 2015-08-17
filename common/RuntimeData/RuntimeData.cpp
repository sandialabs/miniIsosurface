/*
 * RuntimeData.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "RuntimeData.h"

template<typename T>
RuntimeData<T>::RuntimeData() : doc("Marching Cubes", "0.1", ".", "yaml_out.yaml") {
	edgeIndices=0;
}

template<typename T>
RuntimeData<T>::~RuntimeData() {
	delete [] edgeIndices;
}

template<typename T>
void RuntimeData<T>::initEdgeIndices(void) {
	edgeIndices = new EdgeIndexer<T>(ext);
	unsigned mapSize = edgeIndices->nAllEdges / 8; // Very approximate hack..
	pointMap.rehash(mapSize);
}


// Must instantiate class for separate compilation
template class RuntimeData<float_t> ;
