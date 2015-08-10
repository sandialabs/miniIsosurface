/*
 * marchBlock.h
 *
 *  Created on: Jul 10, 2015
 *      Author: sjmunn
 */

#ifndef ALGORITHM_MARCHBLOCK_H_
#define ALGORITHM_MARCHBLOCK_H_

// Common utility headers -------------
// Standard C/C++ library
#include"../../utils/includes.h"

// Data Objects
#include"../../utils/types.h"

// Constants
#include"../../utils/Constants/MarchingCubesTables.h"

// Local implementation headers -------
// Algorithm objects/functions
#include"gradients.h"


static inline float_t lerp(float_t a, float_t b, float_t w);

void extractIsosurfaceFromBlock(const Image3D_t &, const unsigned ext[6],
		float_t, PointMap_t &,EdgeIndexer_t &, TriangleMesh_t *);

#endif /* ALGORITHM_MARCHBLOCK_H_ */
