/*
 * GeneralContext.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "GeneralContext.h"

template<typename T>
GeneralContext<T>::GeneralContext(MarchAlgorithm<T> * implementation) : doc("Marching Cubes", "0.1", ".", "yaml_out.yaml") {
	_strategy=implementation;
}

template<typename T>
GeneralContext<T>::GeneralContext() : doc("Marching Cubes", "0.1", ".", "yaml_out.yaml") {
	_strategy=0;
}

template<typename T>
GeneralContext<T>::~GeneralContext() {
}

template<typename T>
void GeneralContext<T>::setAlgorithm(MarchAlgorithm<T> * implementation) {
	_strategy=implementation;
}

template<typename T>
void GeneralContext<T>::march(void) {
	_strategy->march(*this);
}

// Must instantiate class for separate compilation
template class GeneralContext<float_t>;
