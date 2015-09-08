/*
 * MergeMPAlgo.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "MergeMPAlgo.h"

template<typename T>
MergeMPAlgo<T>::MergeMPAlgo() {
	grainDim=256;
}

template<typename T>
MergeMPAlgo<T>::MergeMPAlgo(unsigned grain) {
	grainDim=grain;
}

template<typename T>
MergeMPAlgo<T>::~MergeMPAlgo() {
	// TODO Auto-generated destructor stub
}

template<typename T>
unsigned MergeMPAlgo<T>::numBlocks(const Range oneDRange) {
	unsigned numberOfBlocks = (oneDRange.end() - oneDRange.begin() + oneDRange.grain() - 1)
			/ oneDRange.grain();
	return numberOfBlocks;
}

template<typename T>
void MergeMPAlgo<T>::march(GeneralContext<T> &data){

	const unsigned *dims = data.imageIn.getDimension();

	CLOG(logYAML) << "Marching cubes algorithm: MergeMP";
	data.doc.add("Marching cubes algorithm", "MergeMP");

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
	this->setGlobalVariables(data);

	#pragma omp parallel
	{
		TriangleMesh_t threadMesh;
		PointMap_t threadPointMap;
		DuplicateRemover threadDuplicateRemover;

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
			//threadPointMap.reserve(mapSize);

			MarchAlgorithm<T>::extractIsosurfaceFromBlock(data.imageIn, blockExtent,
					data.isoval, threadPointMap, *(this->globalEdgeIndices), threadMesh);

			threadDuplicateRemover.setArrays(threadPointMap);
		}

		/*
		 * The next sections merges all the meshes in each block.
		 * It must be done in serial to construct an accurate mesh.
		 */
		#pragma omp critical
		{
			const unsigned nPoints=duplicateRemover.getSize();
			threadDuplicateRemover += nPoints;
			duplicateRemover+=threadDuplicateRemover;
			// The += operator for TriangleMesh3D object is overloaded to merge mesh objects
			meshBeforeMerge += threadMesh;
		}
	}

	duplicateRemover.sortYourSelf();
	duplicateRemover.getNewIndices();

	CLOG(logDEBUG) << "final";
	buildMesh(data.mesh,meshBeforeMerge,duplicateRemover);

	CLOG(logDEBUG1) << "Mesh verts: " << data.mesh.numberOfVertices();
	CLOG(logDEBUG1) << "Mesh tris: " << data.mesh.numberOfTriangles();
}


// Must instantiate class for separate compilation
template class MergeMPAlgo<float_t> ;
