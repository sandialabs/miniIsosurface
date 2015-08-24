/*
 * MpiAlgo.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "MpiAlgo.h"

template<typename T>
MpiAlgo<T>::MpiAlgo(LoadImage3DMPI<T> & inFileHeader, int inPid, int inProcesses) : fileHeader(inFileHeader) {
	pID = inPid;
	processes=inProcesses;
	unsigned maxDim=inFileHeader.getMaxVoumeDimension();

	// Should be
	//grainDim=maxDim/inProcesses;
	grainDim=maxDim/2;
	CLOG(logDEBUG) << "grainDim: " << grainDim;
}

template<typename T>
MpiAlgo<T>::MpiAlgo(unsigned grain) {
	grainDim=grain;
	meshBeforeMerge=0;
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
void MpiAlgo<T>::march(GeneralContext<T> &data) {

//	/*
//	 * MPI stuff
//	 */
//	// Get the number of processes
//	int world_size;
//	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//
//	// Get the rank of the process
//	int world_rank;
//	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//
//	CLOG(logDEBUG) << "Process rank is " << world_rank;

	// DEBUGGER
	{
	    int i = 0;
	    char hostname[256];
	    gethostname(hostname, sizeof(hostname));
	    printf("PID %d on %s ready for attach\n", getpid(), hostname);
	    fflush(stdout);
	    while (0 == i)
	        sleep(5);
	}

	LoadImage3DMPI<float_t> fileData(fileHeader); // We will need multiple data loaders in MPI

	const unsigned *dims = fileHeader.getVolumeDimensions();

	CLOG(logYAML) << "Marching cubes algorithm: MPI";
	data.doc.add("Marching cubes algorithm", "MPI");

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

	TriangleMesh_t processMesh;
	PointMap_t processPointMap;
	DuplicateRemover * processDuplicateRemover=new DuplicateRemover;

	//CLOG(logDEBUG) << "Iteration " << i;
	unsigned blockPageIdx = pID / nblocksPerPage;
	unsigned blockRowIdx = (pID % nblocksPerPage) / numBlockCols;
	unsigned blockColIdx = (pID % nblocksPerPage) % numBlockCols;

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
	fileData.readBlockData(data.imageIn);
	data.imageIn.setToMPIdataBlock();
	data.imageIn.setMPIorigin(blockExtent[0],blockExtent[2],blockExtent[4]);


	unsigned approxNumberOfEdges = 3*(pto-pfrom)*(rto-rfrom)*(cto-cfrom);

	unsigned mapSize = approxNumberOfEdges / 8 + 6; // Very approximate hack..
	//processPointMap.reserve(mapSize);

	MarchAlgorithm<T>::extractIsosurfaceFromBlock(data.imageIn, blockExtent,
			data.isoval, processPointMap, *(this->globalEdgeIndices), processMesh);

	processDuplicateRemover->setArrays(processPointMap);

	processDuplicateRemover->sortYourSelf();
	processDuplicateRemover->getNewIndices();

	buildMesh(data.mesh,processMesh,*processDuplicateRemover);

	CLOG(logDEBUG1) << "Mesh verts: " << data.mesh.numberOfVertices();
	CLOG(logDEBUG1) << "Mesh tris: " << data.mesh.numberOfTriangles();
}


// Must instantiate class for separate compilation
template class MpiAlgo<float_t> ;
