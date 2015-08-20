/*
 * MarchAlgorithm.h
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#ifndef MARCHALGORITHM_H_
#define MARCHALGORITHM_H_

// External includes
#include"../includes.h"

// Reporting Headers
#include"../Reporting/Log.h"

#include"../Constants/MarchingCubesTables.h"
#include"./gradients.h"
#include"./EdgeIndexer.h"

#include"./BlockMarchFunctor.h"

template<typename U> class GeneralContext;

template<typename T>
class MarchAlgorithm {

	typedef Image3D<T> Image3D_type;
	typedef TriangleMesh<T> TriangleMesh_type;
	typedef EdgeIndexer<T> EdgeIndexer_type;
	typedef std::unordered_map<unsigned,unsigned> PointMap_type;

	typedef Triplet<float_t> PositionVector_type;
	typedef Triplet<unsigned> IndexTriplet_type;
public:
	MarchAlgorithm();
	virtual ~MarchAlgorithm();

	virtual void march(GeneralContext<T>&) = 0;

	void setGlobalVariables(GeneralContext<T> &);

	void extractIsosurfaceFromBlock(Image3D_type &, const unsigned [6],
			T, PointMap_type &, EdgeIndexer_type &,
			TriangleMesh_type &);

	EdgeIndexer_type *globalEdgeIndices;
	PointMap_type globalPointMap;
};

#endif /* MARCHALGORITHM_H_ */
