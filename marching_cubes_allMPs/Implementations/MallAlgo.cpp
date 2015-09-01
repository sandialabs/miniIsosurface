/*
 * MallAlgo.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "MallAlgo.h"

template<typename T>
MallAlgo<T>::MallAlgo(LoadImage3DMPI<T> & inFileHeader, unsigned inPid, unsigned inProcesses, Timer * inProcessTimer) : MarchAlgorithm<T>(),fileHeader(inFileHeader) {
	pID = inPid;
	processes=inProcesses;
	const unsigned *dims = fileHeader.getVolumeDimensions();
	grainDimX=dims[0]-1;
	grainDimY=dims[1]-1;

	grainDimZ=(dims[2]-1)/(4*inProcesses);

	if (grainDimZ==0) grainDimZ=1;

	processTimer = inProcessTimer;
}

template<typename T>
MallAlgo<T>::~MallAlgo() {
	// TODO Auto-generated destructor stub
}

template<typename T>
unsigned MallAlgo<T>::numBlocks(const Range oneDRange) {
	unsigned numberOfBlocks = (oneDRange.end() - oneDRange.begin() + oneDRange.grain() - 1)
			/ oneDRange.grain();
	return numberOfBlocks;
}

template<typename T>
bool MallAlgo<T>::testZeroExtent(unsigned * extent) {
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
void MallAlgo<T>::march(GeneralContext<T> &data) {
	const unsigned *dims = fileHeader.getVolumeDimensions();

	CLOG(logWARNING) << "The MPI + OpenMP version of this code is unstable";
	CLOG(logYAML) << "Marching cubes algorithm: MPI + OpenMP";

	data.doc.add("Marching cubes algorithm", "MPI + OpenMP");

	/*
	 * fullRange is the entire image volume
	 */
	Range3D fullRange(0, dims[2] - 1, grainDimZ, 0, dims[1] - 1, grainDimY, 0,
			dims[0] - 1, grainDimX);
	fullRange.extent(data.ext);

	unsigned numBlockPages = this->numBlocks(fullRange.pages());
	unsigned numBlockRows = this->numBlocks(fullRange.rows());
	unsigned numBlockCols = this->numBlocks(fullRange.cols());

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

	unsigned processExtent[6]={dims[0]-1,0,dims[1]-1,0,dims[2]-1,0};

	// Load data for a process
	for (unsigned blockNum=startBlockNum;blockNum<endBlockNum;++blockNum) {
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
//		CLOG(logDEBUG) << "Current Block Extent: " << blockExtent[0] << "  " << blockExtent[1]
//						<< " " << blockExtent[2] << " " << blockExtent[3] << " " << blockExtent[4]
//						<< " " << blockExtent[5];
//		CLOG(logDEBUG1) << "Process extent: " << processExtent[0] << "  " << processExtent[1]
//								<< " " << processExtent[2] << " " << processExtent[3] << " "
//								<< processExtent[4] << " " << processExtent[5];
		for (int i=0;i<5;i+=2) {
			if (blockExtent[i] < processExtent[i]) processExtent[i] = blockExtent[i];
		}
		for (int i=1;i<6;i+=2) {
			if (blockExtent[i] > processExtent[i]) processExtent[i] = blockExtent[i];
		}
	}
	LoadImage3DMPI<float_t> fileData(fileHeader);
//	CLOG(logDEBUG1) << "Process extent: " << processExtent[0] << "  " << processExtent[1]
//						<< " " << processExtent[2] << " " << processExtent[3] << " "
//						<< processExtent[4] << " " << processExtent[5];

	fileData.setBlockExtent(processExtent);
	fileData.readBlockData(data.imageIn);
	data.imageIn.setToMPIdataBlock();
	data.imageIn.setMPIorigin(processExtent[0],processExtent[2],processExtent[4]);

	#pragma omp parallel
	{
		TriangleMesh_t threadMesh;
		PointMap_t threadPointMap;
		DuplicateRemover threadDuplicateRemover;
		unsigned blockId=0;
		EdgeIndexer<T> localEdges(processExtent);

		#pragma omp for nowait
		for (unsigned blockNum=startBlockNum;blockNum<endBlockNum;++blockNum) {
			//CLOG(logINFO) << "Process number " << pID << " is working on " << blockNum << " blocks of " << nblocks;

			//CLOG(logDEBUG) << "Block num " << blockNum;
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

			//unsigned approxNumberOfEdges = 3*(pto-pfrom)*(rto-rfrom)*(cto-cfrom);

			//unsigned mapSize = approxNumberOfEdges / 8 + 6; // Approx # of edges in map

			MarchAlgorithm<T>::extractIsosurfaceFromBlock(data.imageIn, blockExtent,
					data.isoval, threadPointMap, localEdges, threadMesh);

			threadDuplicateRemover.setArrays(threadPointMap);
			blockId=blockNum;
		}

		#pragma omp critical
		{
//			std::string outFile = "threadMesh.";
//			outFile = outFile + std::to_string(static_cast<long long int>(blockId));
//			saveTriangleMesh(&threadMesh, outFile.c_str());
			const unsigned nPoints=processDuplicateRemover.getSize();
			threadDuplicateRemover += nPoints;
			processDuplicateRemover+=threadDuplicateRemover;

			CLOG(logDEBUG1) << "Thread mesh verts: " << threadMesh.numberOfVertices();
			CLOG(logDEBUG1) << "Thread mesh tris: " << threadMesh.numberOfTriangles();
			// The += operator for TriangleMesh3D object is overloaded to merge mesh objects
			processMesh += threadMesh;
		}
	}

	CLOG(logDEBUG1) << "Process mesh verts: " << processMesh.numberOfVertices();
	CLOG(logDEBUG1) << "Process mesh tris: " << processMesh.numberOfTriangles();

	if (processDuplicateRemover.getSize() > 3) {
		processDuplicateRemover.sortYourSelf();
		processDuplicateRemover.getNewIndices();

		buildMesh(data.mesh,processMesh,processDuplicateRemover);
	}

	CLOG(logDEBUG1) << "Mesh verts: " << data.mesh.numberOfVertices();
	CLOG(logDEBUG1) << "Mesh tris: " << data.mesh.numberOfTriangles();
}


// Must instantiate class for separate compilation
template class MallAlgo<float_t> ;
