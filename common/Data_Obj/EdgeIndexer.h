/*
 * TriMeshWithEdges.h
 *
 *  Created on: Jul 28, 2015
 *      Author: sjmunn
 */

#ifndef DATA_OBJ_EDGEINDEXER_H_
#define DATA_OBJ_EDGEINDEXER_H_

#include"../includes.h"

#include "TriangleMesh.h"



// Reporting Headers
#include"../Reporting/Log.h"
#include"../Reporting/RunTime_errors.h"

typedef Triplet<unsigned> Point3dIdx;

template<typename T>
class EdgeIndexer {

	typedef Triplet<T> Point3d;
	typedef Triplet<T> Vector3d;
	typedef Triplet<unsigned> Triangle;
public:
	EdgeIndexer(unsigned *);
	virtual ~EdgeIndexer() {};

	// Adding to points to the edge index
	const unsigned getEdgeIndex(unsigned, unsigned, unsigned, int) const;
	const void getPointCoordinates(unsigned,Point3dIdx&,Point3dIdx&) const;

	// This function is not a member b/c we keep I/O separate
	template <typename Anything> friend void saveTriangleMesh(const EdgeIndexer<Anything> *, const char *);

	unsigned nAllEdges;
	unsigned rangeX, rangeY, rangeZ;
	unsigned nXedges, nXYedges;
private:
	// Set of 3 edge index calculators
	const unsigned edgeIndexXaxis(unsigned, unsigned, unsigned, int) const;
	const unsigned edgeIndexYaxis(unsigned, unsigned, unsigned, int) const;
	const unsigned edgeIndexZaxis(unsigned, unsigned, unsigned, int) const;
};

#endif /* DATA_OBJ_EDGEINDEXER_H_ */
