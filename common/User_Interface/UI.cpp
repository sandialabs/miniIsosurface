/*
 * UI.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: sjmunn
 */

#include "UI.h"

template<typename T>
UI<T>::UI(int argc, char* argv[]) {

	filePath = argv[1];
	outFilePath = argv[2];
	isoval = atof(argv[3]);
	this->checkDebugLevel();
}

template<typename T>
const char * UI<T>::getFile(void) const {
	return filePath;
}

template<typename T>
const char * UI<T>::outFile(void) const {
	return outFilePath;
}

template<typename T>
T UI<T>::getIsoVal(void) const {
	return isoval;
}


template<typename T>
void UI<T>::checkDebugLevel(void) const {
	if (LOG::ReportingLevel() == logDEBUG_Step) {
		CLOG(logWARNING) << "Step by step debugging output is active";
		CLOG(logWARNING)
				<< "Only image volumes with up to 1000 points will be accepted";
	}
}

template<typename T>
void UI<T>::checkArgs(int argc) {
	if (argc != 4) {
		throw wrong_arguments(argc);
	}
}

// Must instantiate class for separate compilation
template class UI<float_t> ;
