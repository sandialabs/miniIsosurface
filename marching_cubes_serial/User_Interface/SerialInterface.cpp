/*
 * MergeMPterface.cpp
 *
 *  Created on: Aug 3, 2015
 *      Author: sjmunn
 */

#include "SerialInterface.h"

template<typename T>
void SerialInterface<T>::marchImplemtation(SerialData<T> &data) {
	SerialAlgo<T> algorithm;
	algorithm.visit(data);
}

// Must instantiate class for separate compilation
template class SerialInterface<float_t> ;
