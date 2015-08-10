/*
 * pointCount.cpp
 *
 *  Created on: Aug 7, 2015
 *      Author: sjmunn
 */

#include"pointCount.h"

static inline void counter(float_t * vals, float_t &isoval, unsigned &nUniquePts, PointMap_t &ptMapDebug, unsigned &edgeIndex) {
	if (vals[0] >= isoval && vals[1] < isoval) {
		if (!(vals[0] == isoval && vals[1] == isoval)) {
			nUniquePts++;
			if (ptMapDebug.find(edgeIndex) == ptMapDebug.end()) {
				// not found -- this is a new point
				ptMapDebug[edgeIndex] = 9;
			}
		}
	}
	else if (vals[0] < isoval && vals[1] >= isoval) {
		if (!(vals[0] == isoval && vals[1] == isoval)) {
			nUniquePts++;
			if (ptMapDebug.find(edgeIndex) == ptMapDebug.end()) {
				// not found -- this is a new point
				ptMapDebug[edgeIndex] = 9;
			}
		}
	}
}

unsigned countPointsInBlock(const Image3D_t &vol, const unsigned ext[6],
		float_t isoval, EdgeIndexer_t &edgeIndices, PointMap_t &ptMapDebug) {
	/*
	 * This function iterates through all the edges (not cubes)
	 * in a data block and returns the number of unique points
	 * there are.
	 */
	unsigned nEdges=edgeIndices.nAllEdges;
	const unsigned *dims = vol.getDimension();
	const float_t *buffer = vol.getData();
	unsigned sliceSize = dims[0] * dims[1];

	unsigned bufferIdx=ext[0]+ ( ext[2] * dims[0]) + (ext[4] * sliceSize);
	unsigned idx;
	float_t val[2];

	unsigned nUniquePts=0;

//	CLOG(logDEBUG) << "nExedges: " << edgeIndices.nXedges;
//	CLOG(logDEBUG) << "nExedges: " << edgeIndices.nXYedges;
//	CLOG(logDEBUG) << "nExedges: " << edgeIndices.nAllEdges;
	/*
	 * For cache efficiency, create two parallel buffers
	 */
	const float_t *X1buffer = &buffer[bufferIdx];
	const float_t *X2buffer = &buffer[bufferIdx + 1];

	// For edges parallel to the x-axis
	for (unsigned iEdge = 0; iEdge<edgeIndices.nXedges; ++iEdge) {
		idx=iEdge+iEdge/edgeIndices.rangeX;
		// get cell-points values
		val[0] = X1buffer[idx];
		val[1] = X2buffer[idx];
		counter(val,isoval,nUniquePts, ptMapDebug,iEdge);
	}

	X2buffer = &buffer[bufferIdx + dims[0]];

	// For edges parallel to the y-axis
	for (unsigned iEdge = edgeIndices.nXedges; iEdge<edgeIndices.nXYedges; ++iEdge) {
		idx=iEdge-edgeIndices.nXedges
				+(edgeIndices.rangeX+1)*((iEdge-edgeIndices.nXedges)
				/((edgeIndices.rangeX+1)*edgeIndices.rangeY));
		// get cell-points values
		val[0] = X1buffer[idx];
		val[1] = X2buffer[idx];
		counter(val,isoval,nUniquePts, ptMapDebug,iEdge);
	}

	X2buffer = &buffer[bufferIdx + sliceSize];

	// For edges parallel to the z-axis
	for (unsigned iEdge = edgeIndices.nXYedges; iEdge<=edgeIndices.nAllEdges; ++iEdge) {
		idx=iEdge-edgeIndices.nXYedges;
		// get cell-points values
		val[0] = X1buffer[idx];
		val[1] = X2buffer[idx];
		counter(val,isoval,nUniquePts, ptMapDebug,iEdge);
	}

	CLOG(logDEBUG1) << "# of unique points: " << nUniquePts;
}

