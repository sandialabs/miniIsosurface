/*
 * Image3D.cpp
 *
 *  Created on: Jul 6, 2015
 *      Author: sjmunn
 */

#include"Image3D.h"

template<typename T>
Image3D<T>::Image3D() :
		npoints(0), data(0) {

	dim[0] = dim[1] = dim[2] = 0;
	spacing[0] = spacing[1] = spacing[2] = 0.0;
	origin[0] = origin[1] = origin[2] = 0.0;
	isMPIdataBlock=false;
	sliceSize=0;
}

template<typename T>
Image3D<T>::~Image3D() {
	delete[] this->data;
}

template<typename T>
void Image3D<T>::setDimension(unsigned xdim, unsigned ydim, unsigned zdim) {
	this->dim[0] = xdim;
	this->dim[1] = ydim;
	this->dim[2] = zdim;
	sliceSize=dim[0]*dim[1];
}

template<typename T>
void Image3D<T>::setSpacing(T xspc, T yspc, T zspc) {
	this->spacing[0] = xspc;
	this->spacing[1] = yspc;
	this->spacing[2] = zspc;
}

template<typename T>
void Image3D<T>::setOrigin(T x, T y, T z) {
	this->origin[0] = x;
	this->origin[1] = y;
	this->origin[2] = z;
}

template<typename T>
void Image3D<T>::setToMPIdataBlock(void) {
	isMPIdataBlock=true;
}

template<typename T>
void Image3D<T>::allocate() {
	this->npoints = dim[0] * dim[1] * dim[2];
	delete[] this->data;
	this->data = new T[npoints];
}

template<typename T>
T* Image3D<T>::getData() {
	return this->data;
}

template<typename T>
const unsigned* Image3D<T>::getDimension() const {
	return this->dim;
}

template<typename T>
unsigned Image3D<T>::getNumberOfPoints() const {
	return this->npoints;
}

template<typename T>
const T* Image3D<T>::getSpacing() const {
	return this->spacing;
}

template<typename T>
const T* Image3D<T>::getOrigin() const {
	return this->origin;
}

template<typename T>
const T* Image3D<T>::getData() const {
	CLOG(logWARNING) << "Image3D<T>::getData() call";
	CLOG(logWARNING) << "Will be deleting this function soon";
	return this->data;
}

template<typename T>
void Image3D<T>::report(YAML_Doc &doc) const {
	// Console output,
	CLOG(logYAML) << "File x-dimension: " << static_cast<long long int>(dim[0]);
	CLOG(logYAML) << "File y-dimension: " << static_cast<long long int>(dim[1]);
	CLOG(logYAML) << "File z-dimension: " << static_cast<long long int>(dim[2]);

	CLOG(logYAML) << "Number of points in image volume: " << static_cast<long long int>(npoints);

	CLOG(logINFO) << "File x-spacing: " << static_cast<long long int>(spacing[0]);
	CLOG(logINFO) << "File y-spacing: " << static_cast<long long int>(spacing[1]);
	CLOG(logINFO) << "File z-spacing: " << static_cast<long long int>(spacing[2]);

	CLOG(logINFO) << "File x-origin: " << static_cast<long long int>(origin[0]);
	CLOG(logINFO) << "File y-origin: " << static_cast<long long int>(origin[1]);
	CLOG(logINFO) << "File z-origin: " << static_cast<long long int>(origin[2]);

	// YAML output,
	doc.add("File x-dimension",static_cast<long long int>(dim[0]));
	doc.add("File y-dimension",static_cast<long long int>(dim[1]));
	doc.add("File z-dimension",static_cast<long long int>(dim[2]));

	doc.add("Number of points in image volume", static_cast<long long int>(npoints));
}

//template<typename T>
//void Image3D<T>::setImage3DOutputBuffers(const unsigned xIdx, const unsigned yIdx, const unsigned zIdx) {
//	bufferIdx= xIdx + (yIdx * dim[0]) + (zIdx * sliceSize);
//
//	X1buffer = &data[bufferIdx];
//	X2buffer = &data[bufferIdx + dim[0]];
//	X3buffer = &data[bufferIdx + sliceSize];
//	X4buffer = &data[bufferIdx + dim[0] + sliceSize];
//}

//template<typename T>
//void Image3D<T>::getVertexValues(T *vertexVals, unsigned xIdx, unsigned xExtent) {
//
//	vertexVals[0] = X1buffer[xIdx-xExtent];
//	vertexVals[1] = X1buffer[xIdx-xExtent+1];
//
//	vertexVals[2] = X2buffer[xIdx-xExtent+1];
//	vertexVals[3] = X2buffer[xIdx-xExtent];
//
//	vertexVals[4] = X3buffer[xIdx-xExtent];
//	vertexVals[5] = X3buffer[xIdx-xExtent+1];
//
//	vertexVals[6] = X4buffer[xIdx-xExtent+1];
//	vertexVals[7] = X4buffer[xIdx-xExtent];
//}

//template<typename T>
//void Image3D<T>::getValsForGradient(T (& x)[3][2], const unsigned xIdxGlobal, const unsigned yIdxGlobal, const unsigned zIdxGlobal) const {
//	/*
//	 * Gradient computation is only done on a per need basis
//	 * Setting up cache buffers for this would not improve performance enough
//	 */
//	// Assuming bufferIdx is updated with the marching, we just need to add distance along the x-axis
//	unsigned ptIdxOnBuffer=xIdxGlobal + yIdxGlobal * dim[0] + zIdxGlobal * sliceSize;
//	if (xIdxGlobal == 0) {
//		x[0][0] = data[ptIdxOnBuffer + 1];
//		x[0][1] = data[ptIdxOnBuffer];
//	} else if (xIdxGlobal == (dim[0] - 1)) {
//		x[0][0] = data[ptIdxOnBuffer];
//		x[0][1] = data[ptIdxOnBuffer - 1];
//	} else {
//		x[0][0] = data[ptIdxOnBuffer + 1];
//		x[0][1] = data[ptIdxOnBuffer - 1];
//	}
//
//	if (yIdxGlobal == 0) {
//		x[1][0] = data[ptIdxOnBuffer + dim[0]];
//		x[1][1] = data[ptIdxOnBuffer];
//	} else if (yIdxGlobal == (dim[1] - 1)) {
//		x[1][0] = data[ptIdxOnBuffer];
//		x[1][1] = data[ptIdxOnBuffer - dim[0]];
//	} else {
//		x[1][0] = data[ptIdxOnBuffer + dim[0]];
//		x[1][1] = data[ptIdxOnBuffer - dim[0]];
//	}
//
//	if (zIdxGlobal == 0) {
//		x[2][0] = data[ptIdxOnBuffer + sliceSize];
//		x[2][1] = data[ptIdxOnBuffer];
//	} else if (zIdxGlobal == (dim[2] - 1)) {
//		x[2][0] = data[ptIdxOnBuffer];
//		x[2][1] = data[ptIdxOnBuffer - sliceSize];
//	} else {
//		x[2][0] = data[ptIdxOnBuffer + sliceSize];
//		x[2][1] = data[ptIdxOnBuffer - sliceSize];
//	}
//}

#include"../Algorithm/Image3DReader.h"
// Must instantiate class for separate compilation
template class Image3D<float_t> ;
