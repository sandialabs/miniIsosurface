/*
 * SerialAlgo.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "SerialAlgo.h"

template<typename T>
serial<T>::SerialAlgo() {
	// TODO Auto-generated constructor stub

}

template<typename T>
serial<T>::~SerialAlgo() {
	// TODO Auto-generated destructor stub
}

template<typename T>
void serial<T>::visit(MCdata<T> *mcDat){
	std::cout << "Hello World" << std::endl;
	MarchAlgorithm<T>::generalMethod();
}
