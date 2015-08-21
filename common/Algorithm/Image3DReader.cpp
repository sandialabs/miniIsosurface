/*
 * Image3DReader.cpp
 *
 *  Created on: Aug 21, 2015
 *      Author: sjmunn
 */

#include "Image3DReader.h"

template<typename T>
Image3DReader<T>::Image3DReader(Image3D<T> * inImageData) {
	bufferIdx=0;
	X1buffer=X2buffer=X3buffer=X4buffer=0;
	imageData=inImageData;
}

template<typename T>
Image3DReader<T>::~Image3DReader() {
	// TODO Auto-generated destructor stub
}

template<typename T>
void Image3DReader<T>::setImage3DOutputBuffers(const unsigned xIdx, const unsigned yIdx, const unsigned zIdx) {
	bufferIdx= xIdx + (yIdx * imageData->dim[0]) + (zIdx * imageData->sliceSize);

	X1buffer = &imageData->data[bufferIdx];
	X2buffer = &imageData->data[bufferIdx + imageData->dim[0]];
	X3buffer = &imageData->data[bufferIdx + imageData->sliceSize];
	X4buffer = &imageData->data[bufferIdx + imageData->dim[0] + imageData->sliceSize];
}


template<typename T>
void Image3DReader<T>::getVertexValues(T *vertexVals, unsigned xIdx, unsigned xExtent) {

	vertexVals[0] = X1buffer[xIdx-xExtent];
	vertexVals[1] = X1buffer[xIdx-xExtent+1];

	vertexVals[2] = X2buffer[xIdx-xExtent+1];
	vertexVals[3] = X2buffer[xIdx-xExtent];

	vertexVals[4] = X3buffer[xIdx-xExtent];
	vertexVals[5] = X3buffer[xIdx-xExtent+1];

	vertexVals[6] = X4buffer[xIdx-xExtent+1];
	vertexVals[7] = X4buffer[xIdx-xExtent];
}

template<typename T>
void Image3DReader<T>::getValsForGradient(T (& x)[3][2], const unsigned xIdxGlobal, const unsigned yIdxGlobal, const unsigned zIdxGlobal) const {
	/*
	 * Gradient computation is only done on a per need basis
	 * Setting up cache buffers for this would not improve performance enough
	 */
	// Assuming bufferIdx is updated with the marching, we just need to add distance along the x-axis
	unsigned ptIdxOnBuffer=xIdxGlobal + yIdxGlobal * imageData->dim[0] + zIdxGlobal * imageData->sliceSize;
	if (xIdxGlobal == 0) {
		x[0][0] = imageData->data[ptIdxOnBuffer + 1];
		x[0][1] = imageData->data[ptIdxOnBuffer];
	} else if (xIdxGlobal == (imageData->dim[0] - 1)) {
		x[0][0] = imageData->data[ptIdxOnBuffer];
		x[0][1] = imageData->data[ptIdxOnBuffer - 1];
	} else {
		x[0][0] = imageData->data[ptIdxOnBuffer + 1];
		x[0][1] = imageData->data[ptIdxOnBuffer - 1];
	}

	if (yIdxGlobal == 0) {
		x[1][0] = imageData->data[ptIdxOnBuffer + imageData->dim[0]];
		x[1][1] = imageData->data[ptIdxOnBuffer];
	} else if (yIdxGlobal == (imageData->dim[1] - 1)) {
		x[1][0] = imageData->data[ptIdxOnBuffer];
		x[1][1] = imageData->data[ptIdxOnBuffer - imageData->dim[0]];
	} else {
		x[1][0] = imageData->data[ptIdxOnBuffer + imageData->dim[0]];
		x[1][1] = imageData->data[ptIdxOnBuffer - imageData->dim[0]];
	}

	if (zIdxGlobal == 0) {
		x[2][0] = imageData->data[ptIdxOnBuffer + imageData->sliceSize];
		x[2][1] = imageData->data[ptIdxOnBuffer];
	} else if (zIdxGlobal == (imageData->dim[2] - 1)) {
		x[2][0] = imageData->data[ptIdxOnBuffer];
		x[2][1] = imageData->data[ptIdxOnBuffer - imageData->sliceSize];
	} else {
		x[2][0] = imageData->data[ptIdxOnBuffer + imageData->sliceSize];
		x[2][1] = imageData->data[ptIdxOnBuffer - imageData->sliceSize];
	}
}

// Must instantiate class for separate compilation
template class Image3DReader<float_t>;
