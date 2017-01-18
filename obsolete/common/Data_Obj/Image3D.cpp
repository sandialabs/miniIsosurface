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
	MPIorigin[0]=MPIorigin[1]=MPIorigin[2]=0;
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
void Image3D<T>::setMPIorigin(unsigned xorg, unsigned yorg, unsigned zorg) {
	this->MPIorigin[0]=xorg;
	this->MPIorigin[1]=yorg;
	this->MPIorigin[2]=zorg;
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

#include"../Algorithm/Image3DReader.h"
// Must instantiate class for separate compilation
template class Image3D<float_t> ;
