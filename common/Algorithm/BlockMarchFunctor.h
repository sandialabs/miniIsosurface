/*
 * BlockMarchFunctor.h
 *
 *  Created on: Aug 19, 2015
 *      Author: sjmunn
 */

#ifndef BLOCKMARCHFUNCTOR_H_
#define BLOCKMARCHFUNCTOR_H_

// External includes
#include"../includes.h"

// Reporting Headers
#include"../Reporting/Log.h"

#include"../Constants/MarchingCubesTables.h"
#include"./gradients.h"
#include"./EdgeIndexer.h"

#include"./gradients.h"

template<typename T>
class BlockMarchFunctor {
	typedef Image3D<T> Image3D_type;
	typedef TriangleMesh<T> TriangleMesh_type;
	typedef EdgeIndexer<T> EdgeIndexer_type;
	typedef std::unordered_map<unsigned,unsigned> PointMap_type;

	typedef Triplet<T> PositionVector_type;
	typedef Triplet<unsigned> IndexTriplet_type;
public:
	BlockMarchFunctor(Image3D_type &vol, const unsigned blockExt[6],
			T isoval, PointMap_type &pointMap, EdgeIndexer_type &edgeIndices,
			TriangleMesh_type &mesh);
	virtual ~BlockMarchFunctor();
private:
	static T lerp(T a, T b, T w);
	void updateBuffers(void);
	static int findCaseId(T*,T);

	// ===== Block parameters ===================
	const unsigned *dims;
	const T *origin;
	const T *spacing;
	const T *buffer;
	unsigned sliceSize;

	// ===== Iteration position identifiers =====
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

	// ===== Cell processing variable ===========
	int cellCaseId;
};

#endif /* BLOCKMARCHFUNCTOR_H_ */
