/*
 * MarchAlgorithm.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "MarchAlgorithm.h"

template<typename T>
MarchAlgorithm<T>::MarchAlgorithm() {
	xidx=0;
	yidx=0;
	zidx=0;
	X1buffer=0;
	X2buffer=0;
	X3buffer=0;
	X4buffer=0;
	bufferIdx=0;
}


template<typename T>
MarchAlgorithm<T>::~MarchAlgorithm() {
	// TODO Auto-generated destructor stub
}

template<typename T>
T MarchAlgorithm<T>::lerp(T a, T b, T w) {
	//return ((1.0 - w) * a) + (w * b);
	return a + (w * (b - a));
}

//template<typename T>
//void MarchAlgorithm<T>::updateBuffers(const T *inBuffer) {
//	X1buffer = &inBuffer[bufferIdx];
//	X2buffer = &inBuffer[bufferIdx + dims[0]];
//	X3buffer = &inBuffer[bufferIdx + sliceSize];
//	X4buffer = &inBuffer[bufferIdx + dims[0] + sliceSize];
//}

template<typename T>
void MarchAlgorithm<T>::extractIsosurfaceFromBlock(RuntimeData<T> * inData, const unsigned blockExt[6]) {

	static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

	const unsigned *dims = inData->imageIn.getDimension();
	const T *origin = inData->imageIn.getOrigin();
	const T *spacing = inData->imageIn.getSpacing();
	const T *buffer = inData->imageIn.getData();

	CLOG(logDEBUG1) << "Extent: " << blockExt[0] << " " << blockExt[1] << " "
			<< blockExt[2] << " " << blockExt[3] << " " << blockExt[4] << " " << blockExt[5];

	unsigned sliceSize = dims[0] * dims[1];

	unsigned ptIdx = inData->mesh.numberOfVertices();

	// march through each cell
	T zpos = origin[2] + (T(blockExt[4]) * spacing[2]);
	for (zidx = blockExt[4]; zidx <= blockExt[5]; ++zidx, zpos += spacing[2]) {

		T ypos = origin[1] + (T(blockExt[2]) * spacing[1]);
		for (yidx = blockExt[2]; yidx <= blockExt[3];
				++yidx, ypos += spacing[1]) {

			bufferIdx=blockExt[0]+ (yidx * dims[0]) + (zidx * sliceSize);

			/*
			 * 4 buffers are created to improve cache efficiency
			 * this improves run time by about .1 seconds
			 */
			X1buffer = &buffer[bufferIdx];
			X2buffer = &buffer[bufferIdx + dims[0]];
			X3buffer = &buffer[bufferIdx + sliceSize];
			X4buffer = &buffer[bufferIdx + dims[0] + sliceSize];

			T xpos = origin[0] + (T(blockExt[0]) * spacing[0]);
			for (xidx = blockExt[0]; xidx <= blockExt[1]; ++xidx, xpos +=
					spacing[0]) {

				T pos[8][3], grad[8][3];
				T val[8];

				// get cell-points values
				val[0] = X1buffer[xidx];
				val[1] = X1buffer[xidx+1];

				val[2] = X2buffer[xidx+1];
				val[3] = X2buffer[xidx];

				val[4] = X3buffer[xidx];
				val[5] = X3buffer[xidx+1];

				val[6] = X4buffer[xidx+1];
				val[7] = X4buffer[xidx];

				int cellId = 0;
				for (int i = 0; i < 8; ++i) {
					cellId |= (val[i] >= inData->isoval) ? caseMask[i] : 0;
				}

				// no intersections
				if (cellId == 0 || cellId == 255) {
					continue;
				}

				// get physical position and gradient of the points
				pos[0][0] = xpos;
				pos[0][1] = ypos;
				pos[0][2] = zpos;

				pos[1][0] = xpos + spacing[0];
				pos[1][1] = ypos;
				pos[1][2] = zpos;

				pos[2][0] = xpos + spacing[0];
				pos[2][1] = ypos + spacing[1];
				pos[2][2] = zpos;

				pos[3][0] = xpos;
				pos[3][1] = ypos + spacing[1];
				pos[3][2] = zpos;

				pos[4][0] = xpos;
				pos[4][1] = ypos;
				pos[4][2] = zpos + spacing[2];

				pos[5][0] = xpos + spacing[0];
				pos[5][1] = ypos;
				pos[5][2] = zpos + spacing[2];

				pos[6][0] = xpos + spacing[0];
				pos[6][1] = ypos + spacing[1];
				pos[6][2] = zpos + spacing[2];

				pos[7][0] = xpos;
				pos[7][1] = ypos + spacing[1];
				pos[7][2] = zpos + spacing[2];

				// get the triangles to generate
				const int *edges = MarchingCubesTables::getCaseTrianglesEdges(
						cellId);
				for (; *edges != -1; edges += 3) {

					unsigned tri[3];
					for (int i = 0; i < 3; ++i) {
						int v1 =
								MarchingCubesTables::getEdgeVertices(edges[i])[0];
						int v2 =
								MarchingCubesTables::getEdgeVertices(edges[i])[1];
						T w = (inData->isoval - val[v1]) / (val[v2] - val[v1]);

						// interpolate vertex position
						T  newPtCoordinates[3];
						PositionVector_type newPt;

						bool exists = false;
						unsigned pointId;

						unsigned edgeIndex = inData->edgeIndices->getEdgeIndex(xidx, yidx,
								zidx, edges[i]);
						if (inData->pointMap.find(edgeIndex) == inData->pointMap.end()) {
							// not found -- this is a new point
							inData->pointMap[edgeIndex] = ptIdx;
							for (int iAxis = 0; iAxis < 3; iAxis++) {
								newPtCoordinates[iAxis] = lerp(pos[v1][iAxis],
										pos[v2][iAxis], w);
							}
							newPt.setCoordinates(newPtCoordinates);
						} else {
							// found -- we already have this point
							exists=true;
							pointId = inData->pointMap[edgeIndex];
						}

						if (!exists) {
							inData->mesh.addPoint(newPt);

							computeAllGradients(xidx, yidx, zidx, buffer, dims,spacing,grad);

							T norm[3];
							for (int iAxis = 0; iAxis < 3; iAxis++) {
								norm[iAxis] = lerp(grad[v1][iAxis],
										grad[v2][iAxis], w);
							}
							inData->mesh.addNormal(norm);

							tri[i] = ptIdx++;
						} else {
							tri[i] = pointId;
						}
					}

					if (tri[0] == tri[1] || tri[1] == tri[2]
							|| tri[2] == tri[0]) {
					} else {
						inData->mesh.addTriangle(tri);
						CLOG(logDEBUG_Step) << "Num tri " << inData->mesh.numberOfTriangles();
					}
				}
			}
		}
	}
}


// Must instantiate class for separate compilation
template class MarchAlgorithm<float_t> ;
