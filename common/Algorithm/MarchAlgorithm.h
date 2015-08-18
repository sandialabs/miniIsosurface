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

	void setGlobalVariables(GeneralContext<T> &inData);
	void extractIsosurfaceFromBlock(GeneralContext<T> *inData, const unsigned blockExt[6]);
private:
	static T lerp(T a, T b, T w);
private:
	// Algorithm iteration position
	unsigned xidx, yidx, zidx;
	unsigned bufferIdx;
	/*
	 * 4 buffers improve cache efficiency
	 * this speeds up run time by about .1 seconds
	 */
	const T *X1buffer;
	const T *X2buffer;
	const T *X3buffer;
	const T *X4buffer;

	EdgeIndexer_type *globalEdgeIndices;
	PointMap_type globalPointMap;
};

#endif /* MARCHALGORITHM_H_ */
