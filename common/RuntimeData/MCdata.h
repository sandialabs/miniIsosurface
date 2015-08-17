/*
 * MCdata.h
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

/*
 * Visitor design pattern ---
 * Visited:
 * MCdata		  [parent]			general marching cubes data objects
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

#ifndef MCDATA_H_
#define MCDATA_H_

// Data Object definitions
#include"../Data_Obj/Image3D.h"
#include"../Data_Obj/TriangleMesh.h"
#include"../Data_Obj/Triplet.h"
#include"../Data_Obj/MapReverse.h"
#include"../Data_Obj/EdgeIndexer.h"

#include"../Algorithm/MarchAlgorithm.h"

template<typename T>
class MCdata : public RuntimeData<T> {
public:
	MCdata();
	virtual ~MCdata();

	void accept(MarchAlgorithm<T> &ma);

private:
	Image3D<T> imageIn;
	unsigned ext[6];
	T isoval;
	std::unordered_map<unsigned,unsigned> pointMap;
	//EdgeIndexer<T> edgeIndices;
	TriangleMesh<T> mesh;
};

#endif /* MCDATA_H_ */
