/*
 * pointCount.cpp
 *
 *  Created on: Aug 7, 2015
 *      Author: sjmunn
 */

#include"pointCount.h"

static inline int counter(float_t * vals, float_t &isoval, unsigned &edgeIndex) {
	if (vals[0] >= isoval && vals[1] < isoval) {
		return 1;
	}
	else if (vals[0] < isoval && vals[1] >= isoval) {
		return 1;
	}
	else {
		return 0;
	}
}

unsigned countPointsInBlock(const Image3D_t &vol, const unsigned ext[6],
		float_t isoval, EdgeIndexer_t &edgeIndices) {
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

	/*
	 * For cache efficiency, create two parallel buffers
	 */
	const float_t *X1buffer = &buffer[bufferIdx];
	const float_t *X2buffer = &buffer[bufferIdx + 1];
	const float_t *Ybuffer = &buffer[bufferIdx + dims[0]];
	const float_t *Zbuffer = &buffer[bufferIdx + sliceSize];
	unsigned nUniquePts=0;

	unsigned iEdge;

	#pragma omp parallel private(val)
	{

		// For edges parallel to the x-axis
		#pragma omp for nowait reduction(+:nUniquePts)
		for (iEdge = 0; iEdge<edgeIndices.nXedges; ++iEdge) {
			idx=iEdge+iEdge/edgeIndices.rangeX;
			// get cell-points values
			val[0] = X1buffer[idx];
			val[1] = X2buffer[idx];
			nUniquePts+=counter(val,isoval,iEdge);
		}

		// For edges parallel to the y-axis
		#pragma omp for nowait reduction(+:nUniquePts)
		for (iEdge = edgeIndices.nXedges; iEdge<edgeIndices.nXYedges; ++iEdge) {
			idx=iEdge-edgeIndices.nXedges
					+(edgeIndices.rangeX+1)*((iEdge-edgeIndices.nXedges)
					/((edgeIndices.rangeX+1)*edgeIndices.rangeY));
			// get cell-points values
			val[0] = X1buffer[idx];
			val[1] = Ybuffer[idx];
			nUniquePts+=counter(val,isoval,iEdge);
		}

		// For edges parallel to the z-axis
		#pragma omp for nowait reduction(+:nUniquePts)
		for (iEdge = edgeIndices.nXYedges; iEdge<=edgeIndices.nAllEdges; ++iEdge) {
			idx=iEdge-edgeIndices.nXYedges;
			// get cell-points values
			val[0] = X1buffer[idx];
			val[1] = Zbuffer[idx];
			nUniquePts+=counter(val,isoval,iEdge);
		}
	}

	CLOG(logDEBUG1) << "# of unique points: " << nUniquePts;
}

