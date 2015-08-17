/*
 * SerialData.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "SerialData.h"

template<typename T>
SerialData<T>::SerialData() {
	// TODO Auto-generated constructor stub

}

template<typename T>
SerialData<T>::~SerialData() {
	// not necessary
}

template<typename T>
void SerialData<T>::accept(MarchAlgorithm<T> &ma) {
	ma.visit(*this);
}

// Must instantiate class for separate compilation
template class SerialData<float_t> ;
