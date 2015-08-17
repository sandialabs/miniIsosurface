/*
 * MapReverse.cpp
 *
 *  Created on: Aug 12, 2015
 *      Author: sjmunn
 */

#include "MapReverse.h"

MapReverse::MapReverse(void) {
	// not necessary
	edge_sorted=false;
	//oldIdx_sorted=false;
}

void MapReverse::preAllocate(unsigned nPoints) {
	dataArray.resize(nPoints);
}

void MapReverse::setArrays(PointMap_type &pointMap) {
	/*
	 * This object reverses the point map
	 * pointMap: edges pointing to point indices
	 * edgeIdxArray : point indices pointing to edge indices
	 */
	unsigned nEdges = pointMap.size();
	this->preAllocate(nEdges);

	// Start the clock
	Timer RunTime;
	//#pragma omp parallel for
	for (unsigned iEdge=0;iEdge<nEdges;++iEdge) {
		unsigned iPoint = pointMap[iEdge];
		dataArray[iPoint].edgeIdx=iEdge;
		dataArray[iPoint].pointIdx=iPoint;
	}

	// Stop Clock
	RunTime.stop();

	//Report and save YAML file
	RunTime.reportTime();
//	for (PointMap_type::const_iterator it=pointMap.begin(); it != pointMap.end(); ++it) {
//		dataArray[it->second].edgeIdx=it->first;
//		dataArray[it->second].pointIdx=it->second;
//	}
}

void MapReverse::sortYourSelf(void) {

	std::sort(dataArray.begin(), dataArray.end(), ByEdgeIdx());
	edge_sorted=true;
	//oldIdx_sorted=false;
}

std::vector<unsigned> MapReverse::oldToNewIdxMap(void) {

	unsigned nPoints=getSize();
	std::vector<unsigned> oldToNewMap;
	oldToNewMap.resize(nPoints);

	#pragma omp parallel for
	for(unsigned iPoint=0;iPoint<nPoints;++iPoint) {
		unsigned ptIdxOld=dataArray[iPoint].pointIdx;
		oldToNewMap[ptIdxOld] = dataArray[iPoint].newPointIdx;
	}

	return oldToNewMap;
}

void MapReverse::getNewIndices(void) {
	if (!edge_sorted) {
		throw object_not_sorted("MapReverse");
	}
	/*
	 * Going through the dataArray and creating new point indices for the data
	 */
	// Initialize run
	unsigned iEdgeIndxPrevious=dataArray[0].edgeIdx;
	unsigned iEdgeIndxCurrent=dataArray[1].edgeIdx;
	dataArray[0].newPointIdx=0;
	unsigned currentPointNum=0;

	unsigned nPoints=getSize();
	for (unsigned iPoint=1; iPoint <= nPoints; iPoint++) {
		iEdgeIndxCurrent=dataArray[iPoint].edgeIdx;
		if (iEdgeIndxCurrent!=iEdgeIndxPrevious) currentPointNum++;

		dataArray[iPoint].newPointIdx=currentPointNum;
		iEdgeIndxPrevious=iEdgeIndxCurrent;
	}
}

MapReverse& MapReverse::operator+=(const unsigned increment) {

	# pragma omp parallel for
	for (unsigned iPoint=0; iPoint<this->getSize(); iPoint++) {
		this->dataArray[iPoint].pointIdx+=increment;
	}

	return *this;
}

MapReverse& MapReverse::operator+=(const MapReverse&threadMap) {

	dataArray.insert(dataArray.end(), threadMap.dataArray.begin(), threadMap.dataArray.end());

	return *this;
}

const unsigned MapReverse::getSize(void) const {

	return dataArray.size();
}

MapReverse::~MapReverse() {
	// TODO Auto-generated destructor stub
}

