/*
 * gradients.cpp
 *
 *  Created on: Aug 6, 2015
 *      Author: sjmunn
 */

#include"gradients.h"

void computeGradient(unsigned xidx, unsigned yidx, unsigned zidx,
		const float_t *buffer, const unsigned dims[3], const float_t spacing[3],
		float_t grad[3]) {
	unsigned xysize = dims[0] * dims[1];
	unsigned ptidx = xidx + yidx * dims[0] + zidx * xysize;

	if (xidx == 0) {
		float_t x1 = buffer[ptidx + 1];
		float_t x2 = buffer[ptidx];
		grad[0] = (x2 - x1) / spacing[0];
	} else if (xidx == (dims[0] - 1)) {
		float_t x1 = buffer[ptidx];
		float_t x2 = buffer[ptidx - 1];
		grad[0] = (x2 - x1) / spacing[0];
	} else {
		float_t x1 = buffer[ptidx + 1];
		float_t x2 = buffer[ptidx - 1];
		grad[0] = (0.5 * (x2 - x1)) / spacing[0];
	}

	if (yidx == 0) {
		float_t y1 = buffer[ptidx + dims[0]];
		float_t y2 = buffer[ptidx];
		grad[1] = (y2 - y1) / spacing[1];
	} else if (yidx == (dims[1] - 1)) {
		float_t y1 = buffer[ptidx];
		float_t y2 = buffer[ptidx - dims[0]];
		grad[1] = (y2 - y1) / spacing[1];
	} else {
		float_t y1 = buffer[ptidx + dims[0]];
		float_t y2 = buffer[ptidx - dims[0]];
		grad[1] = (0.5 * (y2 - y1)) / spacing[1];
	}

	if (zidx == 0) {
		float_t z1 = buffer[ptidx + xysize];
		float_t z2 = buffer[ptidx];
		grad[2] = (z2 - z1) / spacing[2];
	} else if (zidx == (dims[2] - 1)) {
		float_t z1 = buffer[ptidx];
		float_t z2 = buffer[ptidx - xysize];
		grad[2] = (z2 - z1) / spacing[2];
	} else {
		float_t z1 = buffer[ptidx + xysize];
		float_t z2 = buffer[ptidx - xysize];
		grad[2] = (0.5 * (z2 - z1)) / spacing[2];
	}
}

void computeAllGradients (unsigned &xidx, unsigned &yidx, unsigned &zidx, const float_t *buffer, const unsigned * dims,
		const float_t spacing[3], float_t (& grad)[8][3]) {
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


