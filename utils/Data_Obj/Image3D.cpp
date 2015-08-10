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

	//doc.add("File x-spacing",static_cast<long long int>(spacing[0]));
	//doc.add("File y-spacing",static_cast<long long int>(spacing[1]));
	//doc.add("File z-spacing",static_cast<long long int>(spacing[2]));

	//doc.add("File x-origin",static_cast<long long int>(origin[0]));
	//doc.add("File y-origin",static_cast<long long int>(origin[1]));
	//doc.add("File z-origin",static_cast<long long int>(origin[2]));
}

// Must instantiate class for separate compilation
template class Image3D<float_t> ;
