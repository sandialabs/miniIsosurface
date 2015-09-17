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
	isMPIimage=inImageData->isMPIdataBlock;
}

template<typename T>
Image3DReader<T>::~Image3DReader() {
	// TODO Auto-generated destructor stub
}

template<typename T>
void Image3DReader<T>::setImage3DOutputBuffers(unsigned xIdx, unsigned yIdx, unsigned zIdx) {
	if (isMPIimage) {
		xIdx=xIdx-imageData->MPIorigin[0];
		yIdx=yIdx-imageData->MPIorigin[1];
		zIdx=zIdx-imageData->MPIorigin[2];
	}
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
void Image3DReader<T>::getValsForGradient(T (& x)[3][2], T (& run)[3], const unsigned xIdxGlobal, const unsigned yIdxGlobal, const unsigned zIdxGlobal) const {
	unsigned xIdx, yIdx, zIdx;

	if (isMPIimage) {
		xIdx=xIdxGlobal-imageData->MPIorigin[0];
		yIdx=yIdxGlobal-imageData->MPIorigin[1];
		zIdx=zIdxGlobal-imageData->MPIorigin[2];
	}
	else {
		xIdx=xIdxGlobal;
		yIdx=yIdxGlobal;
		zIdx=zIdxGlobal;
	}
	/*
	 * Gradient computation is only done on a per need basis
	 * Setting up cache buffers for this would not improve performance enough
	 */
	unsigned ptIdxOnBuffer=xIdx + yIdx * imageData->dim[0] + zIdx * imageData->sliceSize;
	
	const T *spacing = imageData->getSpacing();
	if (xIdx == 0) {
		x[0][0] = imageData->data[ptIdxOnBuffer + 1];
		x[0][1] = imageData->data[ptIdxOnBuffer];
		run[0] = spacing[0];
	} else if (xIdx == (imageData->dim[0] - 1)) {
		x[0][0] = imageData->data[ptIdxOnBuffer];
		x[0][1] = imageData->data[ptIdxOnBuffer - 1];
		run[0] = spacing[0];
	} else {
		x[0][0] = imageData->data[ptIdxOnBuffer + 1];
		x[0][1] = imageData->data[ptIdxOnBuffer - 1];
		run[0] = 2*spacing[0];
	}

	if (yIdx == 0) {
		x[1][0] = imageData->data[ptIdxOnBuffer + imageData->dim[0]];
		x[1][1] = imageData->data[ptIdxOnBuffer];
		run[1] = spacing[1];
	} else if (yIdx == (imageData->dim[1] - 1)) {
		x[1][0] = imageData->data[ptIdxOnBuffer];
		x[1][1] = imageData->data[ptIdxOnBuffer - imageData->dim[0]];
		run[1] = spacing[1];
	} else {
		x[1][0] = imageData->data[ptIdxOnBuffer + imageData->dim[0]];
		x[1][1] = imageData->data[ptIdxOnBuffer - imageData->dim[0]];
		run[1] = 2*spacing[1];
	}

	if (zIdx == 0) {
		x[2][0] = imageData->data[ptIdxOnBuffer + imageData->sliceSize];
		x[2][1] = imageData->data[ptIdxOnBuffer];
		run[2] = spacing[2];
	} else if (zIdx == (imageData->dim[2] - 1)) {
		x[2][0] = imageData->data[ptIdxOnBuffer];
		x[2][1] = imageData->data[ptIdxOnBuffer - imageData->sliceSize];
		run[2] = spacing[2];
	} else {
		x[2][0] = imageData->data[ptIdxOnBuffer + imageData->sliceSize];
		x[2][1] = imageData->data[ptIdxOnBuffer - imageData->sliceSize];
		run[2] = 2*spacing[2];
	}
}

// Must instantiate class for separate compilation
template class Image3DReader<float_t>;
