/*
 * MarchAlgorithm.h
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

#ifndef MARCHALGORITHM_H_
#define MARCHALGORITHM_H_

// External includes
#include"../includes.h"

// Reporting Headers
#include"../Reporting/Log.h"

// Local Includes
#include"../RuntimeData/RuntimeData.h"

#include"../Constants/MarchingCubesTables.h"


template<typename U> class SerialData;

template<typename T>
class MarchAlgorithm {
	typedef SerialData<T> SerialData_type;

	typedef Image3D<T> Image3D_type;
	typedef TriangleMesh<T> TriangleMesh_type;
	typedef EdgeIndexer<T> EdgeIndexer_type;
	typedef std::unordered_map<unsigned,unsigned> PointMap_type;

	typedef Triplet<float_t> PositionVector_type;
	typedef Triplet<unsigned> IndexTriplet_type;
public:
	MarchAlgorithm();
	virtual ~MarchAlgorithm();

	virtual void visit(SerialData<T> &data) = 0;

	// General Marching Cubes methods
	void computeGradient(unsigned xidx, unsigned yidx, unsigned zidx,
			const T *buffer, const unsigned dims[3], const T spacing[3],
			T grad[3]);
	void computeAllGradients (unsigned &xidx, unsigned &yidx, unsigned &zidx, const T *buffer, const unsigned * dims,
			const T spacing[3], T (& grad)[8][3]);

	static T lerp(T a, T b, T w);

	void extractIsosurfaceFromBlock(const Image3D_type &vol, const unsigned ext[6],
			T isoval, PointMap_type &pointMap, EdgeIndexer_type *edgeIndices,
			TriangleMesh_type *mesh);
};

#endif /* MARCHALGORITHM_H_ */
