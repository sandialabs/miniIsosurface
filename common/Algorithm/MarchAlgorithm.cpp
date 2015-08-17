/*
 * MarchAlgorithm.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "MarchAlgorithm.h"

template<typename T>
MarchAlgorithm<T>::MarchAlgorithm() {
	// TODO Auto-generated constructor stub

}


template<typename T>
MarchAlgorithm<T>::~MarchAlgorithm() {
	// TODO Auto-generated destructor stub
}

// General Marching Cubes members,
template<typename T>
void MarchAlgorithm<T>::computeGradient(unsigned xidx, unsigned yidx, unsigned zidx,
		const T *buffer, const unsigned dims[3], const T spacing[3],
		T grad[3]) {
	unsigned xysize = dims[0] * dims[1];
	unsigned ptidx = xidx + yidx * dims[0] + zidx * xysize;

	if (xidx == 0) {
		T x1 = buffer[ptidx + 1];
		T x2 = buffer[ptidx];
		grad[0] = (x2 - x1) / spacing[0];
	} else if (xidx == (dims[0] - 1)) {
		T x1 = buffer[ptidx];
		T x2 = buffer[ptidx - 1];
		grad[0] = (x2 - x1) / spacing[0];
	} else {
		T x1 = buffer[ptidx + 1];
		T x2 = buffer[ptidx - 1];
		grad[0] = (0.5 * (x2 - x1)) / spacing[0];
	}

	if (yidx == 0) {
		T y1 = buffer[ptidx + dims[0]];
		T y2 = buffer[ptidx];
		grad[1] = (y2 - y1) / spacing[1];
	} else if (yidx == (dims[1] - 1)) {
		T y1 = buffer[ptidx];
		T y2 = buffer[ptidx - dims[0]];
		grad[1] = (y2 - y1) / spacing[1];
	} else {
		T y1 = buffer[ptidx + dims[0]];
		T y2 = buffer[ptidx - dims[0]];
		grad[1] = (0.5 * (y2 - y1)) / spacing[1];
	}

	if (zidx == 0) {
		T z1 = buffer[ptidx + xysize];
		T z2 = buffer[ptidx];
		grad[2] = (z2 - z1) / spacing[2];
	} else if (zidx == (dims[2] - 1)) {
		T z1 = buffer[ptidx];
		T z2 = buffer[ptidx - xysize];
		grad[2] = (z2 - z1) / spacing[2];
	} else {
		T z1 = buffer[ptidx + xysize];
		T z2 = buffer[ptidx - xysize];
		grad[2] = (0.5 * (z2 - z1)) / spacing[2];
	}
}

template<typename T>
void MarchAlgorithm<T>::computeAllGradients (unsigned &xidx, unsigned &yidx, unsigned &zidx, const T *buffer, const unsigned * dims,
		const T spacing[3], T (& grad)[8][3]) {
	computeGradient(xidx, yidx, zidx, buffer, dims,
			spacing, grad[0]);
	computeGradient(xidx + 1, yidx, zidx, buffer, dims,
			spacing, grad[1]);
	computeGradient(xidx + 1, yidx + 1, zidx, buffer,
			dims, spacing, grad[2]);
	computeGradient(xidx, yidx + 1, zidx, buffer, dims,
			spacing, grad[3]);
	computeGradient(xidx, yidx, zidx + 1, buffer, dims,
			spacing, grad[4]);
	computeGradient(xidx + 1, yidx, zidx + 1, buffer,
			dims, spacing, grad[5]);
	computeGradient(xidx + 1, yidx + 1, zidx + 1,
			buffer, dims, spacing, grad[6]);
	computeGradient(xidx, yidx + 1, zidx + 1, buffer,
			dims, spacing, grad[7]);
}
