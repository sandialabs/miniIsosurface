/*
 * RuntimeData.h
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

/*
 * Visitor design pattern ---
 * Visited:
 * RuntimeData	  [parent]			general marching cubes data objects
 * SerialData     [child]			data objects specific to the serial implementation
 * OpenMPData	  [child]			data objects specific to the openMP implementation
 * MergeMPData	  [child]			data objects specific to the openMP with mesh merging implementation
 *
 * Visitors:
 * MarchAlgorithm [parent]			general marching cubes methods
 * SerialAlgo	  [child]			marching cubes methods specific to serial implementation
 * OpenMPAlgo	  [child]			marching cubes methods specific to openMP implementation
 * MergeMPAlgo	  [child]			marching cubes methods specific to openMP with mesh merging implementation
 */

#ifndef RUNTIMEDATA_H_
#define RUNTIMEDATA_H_

#include"../includes.h"

// Data Object definitions
#include"../Data_Obj/Image3D.h"
#include"../Data_Obj/TriangleMesh.h"
#include"../Data_Obj/Triplet.h"
#include"../Data_Obj/MapReverse.h"
#include"../Data_Obj/EdgeIndexer.h"

template<typename U> class MarchAlgorithm;

template<typename T>
class RuntimeData {
public:
	RuntimeData();
	virtual ~RuntimeData();

	// Visitor Design Pattern
	virtual void accept(MarchAlgorithm<T> &iVisitor) = 0;

private:
	// General/minimal marching cubes runtime data
	Image3D<T> imageIn;
	unsigned ext[6];
	T isoval;
	std::unordered_map<unsigned,unsigned> pointMap;
	EdgeIndexer<T> edgeIndices;
	TriangleMesh<T> mesh;
};

#endif /* RUNTIMEDATA_H_ */
