/*
 * marchBlock.cpp
 *
 *  Created on: Jul 10, 2015
 *      Author: sjmunn
 */

#include"marchBlock.h"

static inline float_t lerp(float_t a, float_t b, float_t w) {
	//return ((1.0 - w) * a) + (w * b);
	return a + (w * (b - a));
}

void extractIsosurfaceFromBlock(const Image3D_t &vol, const unsigned ext[6],
		float_t isoval, PointMap_t &pointMap, EdgeIndexer_t &edgeIndices,
		TriangleMesh_t *mesh) {

	static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

	const unsigned *dims = vol.getDimension();
	const float_t *origin = vol.getOrigin();
	const float_t *spacing = vol.getSpacing();
	const float_t *buffer = vol.getData();

	CLOG(logDEBUG1) << "Extent: " << ext[0] << " " << ext[1] << " "
			<< ext[2] << " " << ext[3] << " " << ext[4] << " " << ext[5];

	unsigned sliceSize = dims[0] * dims[1];

	unsigned ptIdx = 0;

	unsigned bufferIdx;

	// march through each cell
	float_t zpos = origin[2] + (float_t(ext[4]) * spacing[2]);
	for (unsigned zidx = ext[4]; zidx <= ext[5]; ++zidx, zpos += spacing[2]) {

		float_t ypos = origin[1] + (float_t(ext[2]) * spacing[1]);
		for (unsigned yidx = ext[2]; yidx <= ext[3];
				++yidx, ypos += spacing[1]) {

			bufferIdx=ext[0]+ (yidx * dims[0]) + (zidx * sliceSize);

			/*
			 * 4 buffers are created to improve cache efficiency
			 * this improves run time by about .1 seconds
			 */
			const float_t *X1buffer = &buffer[bufferIdx];
			const float_t *X2buffer = &buffer[bufferIdx + dims[0]];
			const float_t *X3buffer = &buffer[bufferIdx + sliceSize];
			const float_t *X4buffer = &buffer[bufferIdx + dims[0] + sliceSize];

			float_t xpos = origin[0] + (float_t(ext[0]) * spacing[0]);
			for (unsigned xidx = ext[0]; xidx <= ext[1]; ++xidx, xpos +=
					spacing[0]) {

				float_t pos[8][3], grad[8][3];
				float_t val[8];

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
					cellId |= (val[i] >= isoval) ? caseMask[i] : 0;
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
						float_t w = (isoval - val[v1]) / (val[v2] - val[v1]);

						// interpolate vertex position
						float_t  newPtCoordinates[3];
						PositionVector_t newPt;

						bool exists = false;
						unsigned pointId;

						unsigned edgeIndex = edgeIndices.getEdgeIndex(xidx, yidx,
								zidx, edges[i]);
						if (pointMap.find(edgeIndex) == pointMap.end()) {
							// not found -- this is a new point
							pointMap[edgeIndex] = ptIdx;
							for (int iAxis = 0; iAxis < 3; iAxis++) {
								newPtCoordinates[iAxis] = lerp(pos[v1][iAxis],
										pos[v2][iAxis], w);
							}
							newPt.setCoordinates(newPtCoordinates);
						} else {
							// found -- we already have this point
							exists=true;
							pointId = pointMap[edgeIndex];
						}

						if (!exists) {
							mesh->addPoint(newPt);

							computeAllGradients(xidx, yidx, zidx, buffer, dims,spacing,grad);

							float_t norm[3];
							for (int iAxis = 0; iAxis < 3; iAxis++) {
								norm[iAxis] = lerp(grad[v1][iAxis],
										grad[v2][iAxis], w);
							}
							mesh->addNormal(norm);

							tri[i] = ptIdx++;
						} else {
							tri[i] = pointId;
						}
					}

					if (tri[0] == tri[1] || tri[1] == tri[2]
							|| tri[2] == tri[0]) {
					} else {
						mesh->addTriangle(tri);
						CLOG(logDEBUG_Step) << "Num tri " << mesh->numberOfTriangles();
					}
				}
			}
		}
	}
}
