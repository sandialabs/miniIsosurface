/*
 * MapReverse.h
 *
 *  Created on: Aug 12, 2015
 *      Author: sjmunn
 */

#ifndef IMPLEMENTATIONS_MAPREVERSE_H_
#define IMPLEMENTATIONS_MAPREVERSE_H_

#include"../includes.h"

// Reporting Headers
#include"../Reporting/Log.h"
#include"../Reporting/IO_errors.h"
#include"../Reporting/Timer.h"

typedef std::unordered_map<unsigned,unsigned> PointMap_type;

struct edgePointPair {
    unsigned edgeIdx;
    unsigned pointIdx;
    unsigned newPointIdx;
};

struct ByEdgeIdx {
    bool operator()(const edgePointPair &left, const edgePointPair &right) {
        return left.edgeIdx < right.edgeIdx;
    }
};

struct ByPointIdx {
    bool operator()(const edgePointPair &left, const edgePointPair &right) {
        return left.pointIdx < right.pointIdx;
    }
};

class MapReverse {
public:
	MapReverse();
	virtual ~MapReverse();

	// write
	void preAllocate(const unsigned);
	void setArrays(PointMap_type &);
	void sortYourSelf(void);
	std::vector<unsigned> oldToNewIdxMap(void);
	MapReverse& operator+=(const unsigned);
	MapReverse& operator+=(const MapReverse&);
	void getNewIndices(void);

	//read
	const unsigned getSize(void) const;

	std::vector<edgePointPair> dataArray;
	bool edge_sorted;
	//bool oldIdx_sorted;
};

#endif /* IMPLEMENTATIONS_MAPREVERSE_H_ */
