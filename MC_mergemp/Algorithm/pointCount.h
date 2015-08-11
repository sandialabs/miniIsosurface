/*
 * pointCount.h
 *
 *  Created on: Aug 7, 2015
 *      Author: sjmunn
 */

#ifndef ALGORITHM_POINTCOUNT_H_
#define ALGORITHM_POINTCOUNT_H_

// Common utility headers -------------
// Standard C/C++ library
#include"../../utils/includes.h"

// Data Objects
#include"../../utils/types.h"

typedef Triplet<unsigned> Point3dIdx;

unsigned countPointsInVolume(const Image3D_t &, const unsigned ext[6], float_t,
		EdgeIndexer_t &);

unsigned countPointsInBlock(const Image3D_t &, const unsigned ext[6], float_t,
		EdgeIndexer_t &);

#endif /* ALGORITHM_POINTCOUNT_H_ */
