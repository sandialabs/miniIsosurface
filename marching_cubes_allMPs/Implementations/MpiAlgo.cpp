/*
 * MpiAlgo.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "MpiAlgo.h"

template<typename T>
MpiAlgo<T>::MpiAlgo(LoadBigImage<T> & inFileHeader, int inPid, int inProcesses, Timer * inProcessTimer) : MarchAlgorithm<T>(),fileHeader(inFileHeader) {
	pID = inPid;
	processes=inProcesses;
	unsigned maxDim=inFileHeader.getMaxVoumeDimension();

	// Should be
//	float cubeRootProc = static_cast<float>(inProcesses);
//	cubeRootProc=cbrt(cubeRootProc);
//	int nCbrtProcesses=static_cast<int>(cubeRootProc);
//	if (nCbrtProcesses<2) nCbrtProcesses=2;
	grainDim=maxDim/inProcesses;

	processTimer = inProcessTimer;
}

template<typename T>
MpiAlgo<T>::~MpiAlgo() {
	// TODO Auto-generated destructor stub
}

template<typename T>
unsigned MpiAlgo<T>::numBlocks(const Range oneDRange) {
	unsigned numberOfBlocks = (oneDRange.end() - oneDRange.begin() + oneDRange.grain() - 1)
			/ oneDRange.grain();
	return numberOfBlocks;
}

template<typename T>
bool MpiAlgo<T>::testZeroExtent(unsigned * extent) {
	/*
	 * True if zero extent
	 * False if extent non-zero
	 */
	unsigned pseudoVolume=0;
	for (int i=0;i<6;++i) {
		pseudoVolume+=extent[i];
	}
	if (pseudoVolume == 0) {
		return true;
	}
	else {
		return false;
	}
}

template<typename T>
void MpiAlgo<T>::march(GeneralContext<T> &data) {

//	// DEBUGGER
//	int i = 0;
//	char hostname[256];
//	gethostname(hostname, sizeof(hostname));
//	printf("PID %d on %s ready for attach\n", getpid(), hostname);
//	fflush(stdout);
//	while (0==i) {
//		sleep(5);
//	}

	const unsigned *dims = fileHeader.getVolumeDimensions();

	CLOG(logWARNING) << "This MPI + OpenMP version of this code is unstable";
	CLOG(logYAML) << "Marching cubes algorithm: MPI + OpenMP";

	data.doc.add("Marching cubes algorithm", "MPI + OpenMP");

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

	// Distributing jobs
	unsigned blocksPerProcess;
	unsigned startBlockNum;
	unsigned endBlockNum;
	if (nblocks>processes) {
		blocksPerProcess=nblocks/processes;

		startBlockNum=pID*blocksPerProcess;

		// The last process picks up the extra blocks..
		if (pID==processes-1) {
			endBlockNum=nblocks;
		}
		else {
			endBlockNum=startBlockNum+blocksPerProcess;
		}
	}
	else {
		if (pID<nblocks) {
			startBlockNum=pID;
			endBlockNum=pID+1;
		}
		else {
			startBlockNum = endBlockNum = 0;
		}
	}

	TriangleMesh_t processMesh;
	PointMap_t processPointMap;
	DuplicateRemover processDuplicateRemover;

	#pragma omp parallel
	{
		TriangleMesh_t threadMesh;
		PointMap_t threadPointMap;
		DuplicateRemover threadDuplicateRemover;

		#pragma omp for nowait
		for (unsigned blockNum=startBlockNum;blockNum<endBlockNum;++blockNum) {
			LoadBigImage<float_t> fileData(fileHeader); // We will need multiple data loaders in MPI

			unsigned blockPageIdx = blockNum / nblocksPerPage;
			unsigned blockRowIdx = (blockNum % nblocksPerPage) / numBlockCols;
			unsigned blockColIdx = (blockNum % nblocksPerPage) % numBlockCols;

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
			fileData.setBlockExtent(blockExtent);

			processTimer->pause();
			fileData.readBlockData(data.imageIn);
			processTimer->resume();

			data.imageIn.setToMPIdataBlock();
			data.imageIn.setMPIorigin(blockExtent[0],blockExtent[2],blockExtent[4]);


			unsigned approxNumberOfEdges = 3*(pto-pfrom)*(rto-rfrom)*(cto-cfrom);

			//unsigned mapSize = approxNumberOfEdges / 8 + 6; // Approx # of edges in map

			MarchAlgorithm<T>::extractIsosurfaceFromBlock(data.imageIn, blockExtent,
					data.isoval, threadPointMap, *(this->globalEdgeIndices), threadMesh);
			threadDuplicateRemover.setArrays(threadPointMap);
		}

		#pragma omp critical
		{
			const unsigned nPoints=processDuplicateRemover.getSize();
			threadDuplicateRemover += nPoints;
			processDuplicateRemover+=threadDuplicateRemover;
			// The += operator for TriangleMesh3D object is overloaded to merge mesh objects
			processMesh += threadMesh;
		}
	}

	processDuplicateRemover.sortYourSelf();
	processDuplicateRemover.getNewIndices();

	CLOG(logDEBUG) << "final";
	buildMesh(data.mesh,processMesh,processDuplicateRemover);

	CLOG(logDEBUG1) << "Mesh verts: " << data.mesh.numberOfVertices();
	CLOG(logDEBUG1) << "Mesh tris: " << data.mesh.numberOfTriangles();
}


// Must instantiate class for separate compilation
template class MpiAlgo<float_t> ;
