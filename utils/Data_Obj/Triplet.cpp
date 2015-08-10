/*
 * Triplet.cpp
 *
 *  Created on: Jul 14, 2015
 *      Author: sjmunn
 */

#include"Triplet.h"

// Constructor/Destructors
template<typename T>
Triplet<T>::Triplet() {
	coordinates[0] = 0;
	coordinates[1] = 0;
	coordinates[2] = 0;
}

template<typename T>
Triplet<T>::Triplet(const T a, const T b, const T c) {
	coordinates[0] = a;
	coordinates[1] = b;
	coordinates[2] = c;
}

template<typename T>
Triplet<T>::Triplet(const T * pcoords) {
	coordinates[0] = pcoords[0];
	coordinates[1] = pcoords[1];
	coordinates[2] = pcoords[2];
}

template<typename T>
Triplet<T>::~Triplet() {
	// Not needed
}

// Assignment
template<typename T>
void Triplet<T>::setCoordinates(const T a, const T b, const T c) {
	coordinates[0] = a;
	coordinates[1] = b;
	coordinates[2] = c;
}

template<typename T>
void Triplet<T>::setCoordinates(const T * pcoords) {
	coordinates[0] = pcoords[0];
	coordinates[1] = pcoords[1];
	coordinates[2] = pcoords[2];
}

template<typename T>
void Triplet<T>::operator=(const T* pt) {
	this->setCoordinates(pt);
}

template<typename T>
void Triplet<T>::operator=(const Triplet& pt) {
	if (this != &pt) {
		this->setCoordinates(pt.coordinates);
	}
}

template<typename T>
Triplet<T>& Triplet<T>::operator+=(const T term) {
	coordinates[0] += term;
	coordinates[1] += term;
	coordinates[2] += term;

	return *this;
}

// Comparison
template<typename T>
bool Triplet<T>::operator==(const Triplet & Tri) const {
	return (Tri.coordinates[0] == this->coordinates[0]
			&& Tri.coordinates[1] == this->coordinates[1]
			&& Tri.coordinates[2] == this->coordinates[2]);
}

template<typename T>
bool Triplet<T>::operator==(const T * coords) const {
	return (coords[0] == this->coordinates[0]
			&& coords[1] == this->coordinates[1]
			&& coords[2] == this->coordinates[2]);
}

template<typename T>
bool Triplet<T>::operator!=(const Triplet & Tri) const {
	return !(Tri.coordinates[0] == this->coordinates[0]
			&& Tri.coordinates[1] == this->coordinates[1]
			&& Tri.coordinates[2] == this->coordinates[2]);
}

template<typename T>
bool Triplet<T>::operator!=(const T * coords) const {
	return !(coords[0] == this->coordinates[0]
			&& coords[1] == this->coordinates[1]
			&& coords[2] == this->coordinates[2]);
}

// Retrieve coordinates
template<typename T>
const T * Triplet<T>::getCoordinates(void) const {
	return coordinates;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Triplet<T>& triple)
{
    os << triple.coordinates[0] << " / " << triple.coordinates[1] << " / "  << triple.coordinates[2];
    return os;
}

// Must instantiate class for separate compilation
template class Triplet<float_t> ;
template class Triplet<unsigned> ;
template std::ostream& operator<<(std::ostream& os, const Triplet<float_t>& triple);
