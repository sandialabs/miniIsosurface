/*
 * openmp.cpp
 *
 *  Created on: Jul 10, 2015
 *      Author: sjmunn
 */

#include"openmp.h"

static const unsigned grainDim = 256;

namespace openmp {

void extractIsosurface(const Image3D_t &vol, float_t isoval,
		TriangleMesh_t *&mesh, YAML_Doc& doc) {

	TriangleMesh_t beforeMergeMesh;

	const unsigned *dims = vol.getDimension();
	const float_t *origin = vol.getOrigin();
	const float_t *spacing = vol.getSpacing();

	CLOG(logYAML) << "Marching cubes algorithm: OpenMP";
	doc.add("Marching cubes algorithm", "OpenMP");

	float_t range[6];
	for (int i = 0; i < 3; ++i) {
		range[i * 2] = origin[i];
		range[(i * 2) + 1] = origin[i]
				+ (static_cast<float_t>(dims[i] - 1) * spacing[i]);
	}

	/*
	 * fullRange is the entire image volume
	 */
	Range3D fullRange(0, dims[2] - 1, grainDim, 0, dims[1] - 1, grainDim, 0,
			dims[0] - 1, grainDim);

	unsigned fullExtent[6];
	fullRange.extent(fullExtent);
	EdgeIndexer_t edgeIndices(fullExtent);
	MapReverse mapReverse;

	unsigned numBlockPages = numBlocks(fullRange.pages());
	unsigned numBlockRows = numBlocks(fullRange.rows());
	unsigned numBlockCols = numBlocks(fullRange.cols());

	unsigned nblocks = numBlockPages * numBlockRows * numBlockCols;
	unsigned nblocksPerPage = numBlockRows * numBlockCols;
	CLOG(logDEBUG1) << "Number of OpenMP Blocks " << nblocks;

	#pragma omp parallel
	{
		TriangleMesh_t threadMesh;
		PointMap_t threadPointMap;
		MapReverse threadMapReverse;

		#pragma omp for nowait
		for (unsigned i = 0; i < nblocks; ++i) {
			unsigned blockPageIdx = i / nblocksPerPage;
			unsigned blockRowIdx = (i % nblocksPerPage) / numBlockCols;
			unsigned blockColIdx = (i % nblocksPerPage) % numBlockCols;

			unsigned pfrom = blockPageIdx * fullRange.pages().grain();
			unsigned pto = std::min(pfrom + fullRange.pages().grain(),
					fullRange.pages().end());
			unsigned rfrom = blockRowIdx * fullRange.rows().grain();
			unsigned rto = std::min(rfrom + fullRange.rows().grain(),
					fullRange.rows().end());
			unsigned cfrom = blockColIdx * fullRange.cols().grain();
			unsigned cto = std::min(cfrom + fullRange.cols().grain(),
					fullRange.cols().end());

			Range3D blockRange(pfrom, pto, rfrom, rto, cfrom, cto);
			unsigned blockExtent[6];
			blockRange.extent(blockExtent);

			EdgeIndexer_t blockEdgeIndices(blockExtent);
			unsigned approxNumberOfEdges = 3*(pto-pfrom)*(rto-rfrom)*(cto-cfrom);

			unsigned mapSize = approxNumberOfEdges / 8 + 6; // Very approximate hack..
			threadPointMap.rehash(mapSize);

			extractIsosurfaceFromBlock(vol, blockExtent, isoval, threadPointMap,edgeIndices,
					&threadMesh);

			//Timer subTime;
			threadMapReverse.setArrays(threadPointMap);
//			// Stop Clock
//			subTime.stop();
//
//			//Report and save YAML file
//			subTime.reportTime();
		}

		/*
		 * The next sections merges all the meshes in each block.
		 * It must be done in serial to construct an accurate mesh.
		 */
		#pragma omp critical
		{
			const unsigned nPoints=mapReverse.getSize();
			threadMapReverse += nPoints;
			mapReverse+=threadMapReverse;
			// The += operator for TriangleMesh3D object is overloaded to merge mesh objects
			beforeMergeMesh += threadMesh;
		}
	}

	mapReverse.sortYourSelf();
	mapReverse.getNewIndices();

	CLOG(logDEBUG) << "new";
	mesh = new TriangleMesh_t();
	mesh->buildMesh(beforeMergeMesh,mapReverse);

	CLOG(logDEBUG1) << "OpenMP raw mesh verts: " << beforeMergeMesh.numberOfVertices();
	CLOG(logDEBUG1) << "Mesh verts: " << mesh->numberOfVertices();
	CLOG(logDEBUG1) << "Mesh tris: " << mesh->numberOfTriangles();
}

unsigned numBlocks(const Range oneDRange) {
	unsigned numberOfBlocks = (oneDRange.end() - oneDRange.begin() + oneDRange.grain() - 1)
			/ oneDRange.grain();
	return numberOfBlocks;
}

}
