/*
 * DuplicateRemover.h
 *
 *  Created on: Aug 12, 2015
 *      Author: sjmunn
 */

#ifndef IMPLEMENTATIONS_MAPREVERSE_H_
#define IMPLEMENTATIONS_MAPREVERSE_H_

#include"../../common/includes.h"

// Reporting Headers
#include"../../common/Reporting/Log.h"
#include"../../common/Reporting/IO_errors.h"
#include"../../common/Reporting/Timer.h"

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

class DuplicateRemover {
public:
	DuplicateRemover();
//	virtual ~DuplicateRemover();

	// write
	void preAllocate(const unsigned);
	void setArrays(const PointMap_type &);
	void sortYourSelf(void);
	std::vector<unsigned> oldToNewIdxMap(void);
	DuplicateRemover& operator+=(const unsigned);
	DuplicateRemover& operator+=(const DuplicateRemover&);
	void getNewIndices(void);

	//read
	const unsigned getSize(void) const;

	std::vector<edgePointPair> dataArray;
	bool edge_sorted;
	//bool oldIdx_sorted;
};

#endif /* IMPLEMENTATIONS_MAPREVERSE_H_ */
