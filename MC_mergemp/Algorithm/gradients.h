/*
 * gradients.h
 *
 *  Created on: Aug 6, 2015
 *      Author: sjmunn
 */

#ifndef ALGORITHM_GRADIENTS_H_
#define ALGORITHM_GRADIENTS_H_

// Common utility headers -------------
// Standard C/C++ library
#include"../../utils/includes.h"

// Data Objects
#include"../../utils/types.h"

// Constants
#include"../../utils/Constants/MarchingCubesTables.h"

static void computeGradient(unsigned xidx, unsigned yidx, unsigned zidx,
		const float_t *buffer, const unsigned dims[3], const float_t spacing[3],
		float_t grad[3]);

void computeAllGradients (unsigned &xidx, unsigned &yidx, unsigned &zidx, const float_t *buffer, const unsigned * dims,
		const float_t spacing[3], float_t (& grad)[8][3]);

#endif /* ALGORITHM_GRADIENTS_H_ */
