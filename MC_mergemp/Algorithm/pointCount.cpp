/*
 * pointCount.cpp
 *
 *  Created on: Aug 7, 2015
 *      Author: sjmunn
 */

#include"pointCount.h"

static inline int counter(float_t * vals, float_t &isoval) {
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

unsigned countPointsInVolume(const Image3D_t &vol, const unsigned ext[6],
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
			nUniquePts+=counter(val,isoval);
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
			nUniquePts+=counter(val,isoval);
		}

		// For edges parallel to the z-axis
		#pragma omp for nowait reduction(+:nUniquePts)
		for (iEdge = edgeIndices.nXYedges; iEdge<=edgeIndices.nAllEdges; ++iEdge) {
			idx=iEdge-edgeIndices.nXYedges;
			// get cell-points values
			val[0] = X1buffer[idx];
			val[1] = Zbuffer[idx];
			nUniquePts+=counter(val,isoval);
		}
	}

	CLOG(logDEBUG1) << "# of unique points: " << nUniquePts;
}

//unsigned countPointsInBlock(const Image3D_t &vol, const unsigned ext[6],
//		float_t isoval, EdgeIndexer_t &edgeIndices) {
//	/*
//	 * This function iterates through all the edges (not cubes)
//	 * in a data block and returns the number of unique points
//	 * there are.
//	 */
//	CLOG(logDEBUG1) << "Extent: " << ext[0] << " " << ext[1] << " "
//			<< ext[2] << " " << ext[3] << " " << ext[4] << " " << ext[5];
//
//	unsigned nEdges=edgeIndices.nAllEdges;
//	const unsigned *dims = vol.getDimension();
//	const float_t *buffer = vol.getData();
//	unsigned sliceSize = dims[0] * dims[1];
//
//	//unsigned bufferIdx=ext[0]+ ( ext[2] * dims[0]) + (ext[4] * sliceSize);
//	unsigned bufferIdx=0;
//	unsigned idx;
//	float_t val[2];
//
//	/*
//	 * For cache efficiency, create two parallel buffers
//	 */
//	const float_t *X1buffer = &buffer[bufferIdx];
//	const float_t *X2buffer = &buffer[bufferIdx + 1];
//	const float_t *Ybuffer = &buffer[bufferIdx + dims[0]];
//	const float_t *Zbuffer = &buffer[bufferIdx + sliceSize];
//	unsigned nUniquePts=0;
//
//	Point3dIdx pt1,pt2;
//	Point3dIdx originPt(ext[0],ext[2],ext[4]);
//
//	unsigned iEdge;
//
////	#pragma omp parallel private(val,pt1,pt2)
////	{
//
//		// For edges parallel to the x-axis
////		#pragma omp for nowait reduction(+:nUniquePts)
//		for (iEdge = 0; iEdge<edgeIndices.nXedges; ++iEdge) {
//			idx=iEdge;
//			unsigned xoffset = ext[0]+ ( ext[2] * dims[0]) + (ext[4] * sliceSize) + (dims[0]-edgeIndices.rangeX)*(iEdge/edgeIndices.rangeX);
//			idx=idx+xoffset;
//			// get cell-points values
//			CLOG(logDEBUG) << "Sterfs";
//			CLOG(logDEBUG) << iEdge/(edgeIndices.rangeX);
//			CLOG(logDEBUG) << dims[0];
//			CLOG(logDEBUG_Step) <<"x1 - idx: " << idx;
//			CLOG(logDEBUG_Step) <<"x2 - idx: " << idx+1;
//			val[0] = X1buffer[idx];
//			val[1] = X2buffer[idx];
//			nUniquePts+=counter(val,isoval);
//		}
//
//		// For edges parallel to the y-axis
////		#pragma omp for nowait reduction(+:nUniquePts)
//		for (iEdge = edgeIndices.nXedges; iEdge<edgeIndices.nXYedges; ++iEdge) {
//			edgeIndices.getPointCoordinates(iEdge,pt1,pt2);
//			pt1+=originPt;
//			idx=pt1.getCoordinates()[0]
//									 + pt1.getCoordinates()[1]*dims[0]
//									 + pt1.getCoordinates()[2]*sliceSize;
//
//			CLOG(logDEBUG_Step) <<"y1 - idx: " << idx;
//			CLOG(logDEBUG_Step) <<"y2 - idx: " << idx+dims[0];
//			//CLOG(logDEBUG_Step) <<"pt1 : " << pt1.getCoordinates()[0] << " " << pt1.getCoordinates()[1] << " " << pt1.getCoordinates()[2];
//			// get cell-points values
//			val[0] = X1buffer[idx];
//			val[1] = Ybuffer[idx];
//			nUniquePts+=counter(val,isoval);
//		}
//
//		// For edges parallel to the z-axis
////		#pragma omp for nowait reduction(+:nUniquePts)
//		for (iEdge = edgeIndices.nXYedges; iEdge<=edgeIndices.nAllEdges; ++iEdge) {
//			edgeIndices.getPointCoordinates(iEdge,pt1,pt2);
//			pt1+=originPt;
//			idx=pt1.getCoordinates()[0]
//									 + pt1.getCoordinates()[1]*dims[0]
//									 + pt1.getCoordinates()[2]*sliceSize;
//
//			CLOG(logDEBUG_Step) <<"z1 - idx: " << idx;
//			CLOG(logDEBUG_Step) <<"z2 - idx: " << idx+sliceSize;
////			CLOG(logDEBUG_Step) <<"pt1 : " << pt1.getCoordinates()[0] << " " << pt1.getCoordinates()[1] << " " << pt1.getCoordinates()[2];
//			// get cell-points values
//			val[0] = X1buffer[idx];
//			val[1] = Zbuffer[idx];
//			nUniquePts+=counter(val,isoval);
//		}
////	}
//
//	CLOG(logINFO) << "# of unique points: " << nUniquePts;
//	return nUniquePts;
//}

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

	//unsigned bufferIdx=ext[0]+ ( ext[2] * dims[0]) + (ext[4] * sliceSize);
	unsigned bufferIdx=0;
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

	Point3dIdx pt1,pt2;
	Point3dIdx originPt(ext[0],ext[2],ext[4]);

	unsigned iEdge;

	#pragma omp parallel private(val,pt1,pt2)
	{

		// For edges parallel to the x-axis
		#pragma omp for nowait reduction(+:nUniquePts)
		for (iEdge = 0; iEdge<edgeIndices.nXedges; ++iEdge) {
			edgeIndices.getPointCoordinates(iEdge,pt1,pt2);
			pt1+=originPt;
			idx=pt1.getCoordinates()[0]
									 + pt1.getCoordinates()[1]*dims[0]
									 + pt1.getCoordinates()[2]*sliceSize;
			// get cell-points values
			CLOG(logDEBUG_Step) <<"x1 - idx: " << idx;
			CLOG(logDEBUG_Step) <<"x2 - idx: " << idx+1;
			val[0] = X1buffer[idx];
			val[1] = X2buffer[idx];
			nUniquePts+=counter(val,isoval);
		}

		// For edges parallel to the y-axis
		#pragma omp for nowait reduction(+:nUniquePts)
		for (iEdge = edgeIndices.nXedges; iEdge<edgeIndices.nXYedges; ++iEdge) {
			edgeIndices.getPointCoordinates(iEdge,pt1,pt2);
			pt1+=originPt;
			idx=pt1.getCoordinates()[0]
									 + pt1.getCoordinates()[1]*dims[0]
									 + pt1.getCoordinates()[2]*sliceSize;

			CLOG(logDEBUG_Step) <<"y1 - idx: " << idx;
			CLOG(logDEBUG_Step) <<"y2 - idx: " << idx+dims[0];
			//CLOG(logDEBUG_Step) <<"pt1 : " << pt1.getCoordinates()[0] << " " << pt1.getCoordinates()[1] << " " << pt1.getCoordinates()[2];
			// get cell-points values
			val[0] = X1buffer[idx];
			val[1] = Ybuffer[idx];
			nUniquePts+=counter(val,isoval);
		}

		// For edges parallel to the z-axis
		#pragma omp for nowait reduction(+:nUniquePts)
		for (iEdge = edgeIndices.nXYedges; iEdge<=edgeIndices.nAllEdges; ++iEdge) {
			edgeIndices.getPointCoordinates(iEdge,pt1,pt2);
			pt1+=originPt;
			idx=pt1.getCoordinates()[0]
									 + pt1.getCoordinates()[1]*dims[0]
									 + pt1.getCoordinates()[2]*sliceSize;

			CLOG(logDEBUG_Step) <<"z1 - idx: " << idx;
			CLOG(logDEBUG_Step) <<"z2 - idx: " << idx+sliceSize;
//			CLOG(logDEBUG_Step) <<"pt1 : " << pt1.getCoordinates()[0] << " " << pt1.getCoordinates()[1] << " " << pt1.getCoordinates()[2];
			// get cell-points values
			val[0] = X1buffer[idx];
			val[1] = Zbuffer[idx];
			nUniquePts+=counter(val,isoval);
		}
	}
	return nUniquePts;
}

