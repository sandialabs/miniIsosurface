/*
 * MCdata.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#include "MCdata.h"

template<typename T>
MCdata<T>::MCdata() {
	// TODO Auto-generated constructor stub

}

template<typename T>
MCdata<T>::~MCdata() {
	// TODO Auto-generated destructor stub
}

template<typename T>
void MCdata<T>::accept(MarchAlgorithm<T> &ma) {
	ma.visit(this);
}
