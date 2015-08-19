/*
 * OpenMPAlgo.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "OpenMPAlgo.h"

template<typename T>
OpenMPAlgo<T>::OpenMPAlgo() {
	// TODO Auto-generated constructor stub
	grainDim=256;
}

template<typename T>
OpenMPAlgo<T>::OpenMPAlgo(unsigned grain) {
	// TODO Auto-generated constructor stub
	grainDim=grain;
}

template<typename T>
OpenMPAlgo<T>::~OpenMPAlgo() {
	// TODO Auto-generated destructor stub
}

template<typename T>
unsigned OpenMPAlgo<T>::numBlocks(const Range oneDRange) {
	unsigned numberOfBlocks = (oneDRange.end() - oneDRange.begin() + oneDRange.grain() - 1)
			/ oneDRange.grain();
	return numberOfBlocks;
}

template<typename T>
void OpenMPAlgo<T>::march(GeneralContext<T> &data){

	const unsigned *dims = data.imageIn.getDimension();

	CLOG(logYAML) << "Marching cubes algorithm: OpenMP";
	data.doc.add("Marching cubes algorithm", "OpenMP");

	unsigned grainDim=256;

	/*
	 * fullRange is the entire image volume
	 */
	Range3D fullRange(0, dims[2] - 1, grainDim, 0, dims[1] - 1, grainDim, 0,
			dims[0] - 1, grainDim);
	fullRange.extent(data.ext);

	unsigned numBlockPages = numBlocks(fullRange.pages());
	unsigned numBlockRows = numBlocks(fullRange.rows());
	unsigned numBlockCols = numBlocks(fullRange.cols());

	unsigned nblocks = numBlockPages * numBlockRows * numBlockCols;
	unsigned nblocksPerPage = numBlockRows * numBlockCols;
	CLOG(logDEBUG1) << "Number of OpenMP Blocks " << nblocks;
	setGlobalVariables(data);

	#pragma omp parallel
	{
		TriangleMesh_t threadMesh;
		PointMap_t threadPointMap;

		#pragma omp for nowait
		for (unsigned i = 0; i < nblocks; ++i) {
			//CLOG(logDEBUG) << "Iteration " << i;
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

			unsigned approxNumberOfEdges = 3*(pto-pfrom)*(rto-rfrom)*(cto-cfrom);

			unsigned mapSize = approxNumberOfEdges / 8 + 6; // Very approximate hack..
			threadPointMap.rehash(mapSize);

			MarchAlgorithm<T>::extractIsosurfaceFromBlock(data.imageIn, blockExtent,
					data.isoval, threadPointMap, *(this->globalEdgeIndices), threadMesh);
		}

		/*
		 * The next sections merges all the meshes in each block.
		 * It must be done in serial to construct an accurate mesh.
		 */
		#pragma omp critical
		{
			CLOG(logDEBUG1) << "thread Mesh verts: " << threadMesh.numberOfVertices();
			CLOG(logDEBUG1) << "thread Mesh tris: " << threadMesh.numberOfTriangles();
			// The += operator for TriangleMesh3D object is overloaded to merge mesh objects
			data.mesh += threadMesh;
		}
	}
	CLOG(logDEBUG1) << "Mesh verts: " << data.mesh.numberOfVertices();
	CLOG(logDEBUG1) << "Mesh tris: " << data.mesh.numberOfTriangles();
}


// Must instantiate class for separate compilation
template class OpenMPAlgo<float_t> ;
