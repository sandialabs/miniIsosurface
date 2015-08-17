/*
 * RuntimeData.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "RuntimeData.h"

template<typename T>
RuntimeData<T>::RuntimeData() {
	// TODO Auto-generated constructor stub
	edgeIndices=0;
}

template<typename T>
RuntimeData<T>::~RuntimeData() {
	// TODO Auto-generated destructor stub
}

template<typename T>
void RuntimeData<T>::initEdgeIndices(void) {
	edgeIndices = new EdgeIndexer<T>(ext);
}


// Must instantiate class for separate compilation
template class RuntimeData<float_t> ;
