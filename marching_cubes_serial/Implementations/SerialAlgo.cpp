/*
 * SerialAlgo.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "SerialAlgo.h"

template<typename T>
SerialAlgo<T>::SerialAlgo() {
	// TODO Auto-generated constructor stub

}

template<typename T>
SerialAlgo<T>::~SerialAlgo() {
	// TODO Auto-generated destructor stub
}

template<typename T>
void SerialAlgo<T>::visit(SerialData<T> &data){
	std::cout << "Hello World" << std::endl;
}


// Must instantiate class for separate compilation
template class SerialAlgo<float_t> ;
