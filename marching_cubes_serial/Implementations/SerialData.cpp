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
	// TODO Auto-generated destructor stub
}

template<typename T>
void SerialData<T>::accept(MarchAlgorithm<T> &ma) {
	ma.visit(this);
}