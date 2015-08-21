/*
 * MarchAlgorithm.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "MarchAlgorithm.h"

template<typename T>
MarchAlgorithm<T>::MarchAlgorithm() {
	globalEdgeIndices=0;
}

template<typename T>
MarchAlgorithm<T>::~MarchAlgorithm() {
	// TODO Auto-generated destructor stub
}

template<typename T>
void MarchAlgorithm<T>::setGlobalVariables(GeneralContext<T> &inData) {
	globalEdgeIndices = new EdgeIndexer<T>(inData.ext);
	unsigned mapSize = globalEdgeIndices->nAllEdges / 8; // Very approximate hack..
	//globalPointMap.rehash(mapSize);
}

template<typename T>
void MarchAlgorithm<T>::extractIsosurfaceFromBlock(Image3D_type &vol, const unsigned blockExt[6],
		T isoval, PointMap_type &pointMap, EdgeIndexer_type &edgeIndices,
		TriangleMesh_type &mesh) {

	Image3DReader<T> volReader(&vol);
	BlockMarchFunctor<T> marchingCubes(volReader, blockExt, isoval, pointMap, edgeIndices, mesh);
}

#include"../GeneralContext/GeneralContext.h"
// Must instantiate class for separate compilation
template class MarchAlgorithm<float_t> ;
