/*
 * serial.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: sjmunn
 */

#include"serial.h"

static const unsigned grainDim = 4;

namespace serial {

void extractIsosurface(const Image3D_t &vol, float_t isoval,
		TriangleMesh_t *&mesh, YAML_Doc& doc) {

	const unsigned *dims = vol.getDimension();
	const float_t *origin = vol.getOrigin();
	const float_t *spacing = vol.getSpacing();

	CLOG(logYAML) << "Marching cubes algorithm: SERIAL";
	doc.add("Marching cubes algorithm", "SERIAL");

	float_t range[6];
	for (int i = 0; i < 3; ++i) {
		range[i * 2] = origin[i];
		range[(i * 2) + 1] = origin[i]
				+ (static_cast<float_t>(dims[i] - 1) * spacing[i]);
	}

	Range3D cellRange(0, dims[2] - 1, grainDim, 0, dims[1] - 1, grainDim, 0,
			dims[0] - 1, grainDim);
	unsigned extent[6];
	cellRange.extent(extent);

	mesh = new TriangleMesh_t();
	PointMap_t pointMap;

	EdgeIndexer_t edgeIndices(extent);

	unsigned mapSize = edgeIndices.nAllEdges / 8; // Very approximate hack..
	pointMap.rehash(mapSize);

	extractIsosurfaceFromBlock(vol, extent, isoval, pointMap, edgeIndices, mesh);
}

}
