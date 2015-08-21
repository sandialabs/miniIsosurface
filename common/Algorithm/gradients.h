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
#include"../../common/includes.h"

// Data Objects
#include"../../common/types.h"

// Constants
#include"../../common/Constants/MarchingCubesTables.h"

template<typename T>
static void computeGradient(unsigned xidx, unsigned yidx, unsigned zidx, const Image3D<T>*,	T grad[3]);

template<typename T>
void computeAllGradients(unsigned &xidx, unsigned &yidx, unsigned &zidx, const Image3D<T>*, T (& grad)[8][3]);

#endif /* ALGORITHM_GRADIENTS_H_ */
