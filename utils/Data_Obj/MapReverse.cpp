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
	oldIdx_sorted=false;
}

void MapReverse::preAllocate(unsigned nPoints) {
	dataArray.resize(nPoints);
}

void MapReverse::setArrays(const PointMap_type &pointMap) {
	/*
	 * This object reverses the point map
	 * pointMap: edges pointing to point indices
	 * edgeIdxArray : point indices pointing to edge indices
	 */
	unsigned nPoints = pointMap.size();
	this->preAllocate(nPoints);
	#pragma omp parallel for nowait
	for (PointMap_type::const_iterator it=pointMap.begin(); it != pointMap.end(); ++it) {
		dataArray[it->second].edgeIdx=it->first;
		dataArray[it->second].pointIdx=it->second;
	}
}

void MapReverse::sortYourSelf(void) {

	std::sort(dataArray.begin(), dataArray.end(), ByEdgeIdx());
	edge_sorted=true;
	oldIdx_sorted=false;
}

void MapReverse::sortByOldIndex(void) {

	std::sort(dataArray.begin(), dataArray.end(), ByPointIdx());
	edge_sorted=false;
	oldIdx_sorted=true;
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

	# pragma omp parallel for nowait
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

