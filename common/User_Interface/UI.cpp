/*
 * UI.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: sjmunn
 */

#include "UI.h"

template<typename T>
UI<T>::UI(int argc, char* argv[]) {

	if (argc==4) {
		filePath = argv[1];
		outFilePath = argv[2];
		isoval = atof(argv[3]);
		this->checkDebugLevel();
		dataFilePath=0;
		headerSeperate=false;
	}
	else if (argc == 5) {
		filePath = argv[1];
		dataFilePath = argv[2];
		outFilePath = argv[3];
		isoval = atof(argv[4]);
		this->checkDebugLevel();
		headerSeperate=true;
	}

}

template<typename T>
const char * UI<T>::getFile(void) const {
	return filePath;
}

template<typename T>
const char * UI<T>::getDataFile(void) const {
	if (headerSeperate) {
		return dataFilePath;
	}
	else {
		return 0;
	}
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
	if (argc != 4 && argc != 5) {
		throw wrong_arguments(argc);
	}
}

// Must instantiate class for separate compilation
template class UI<float_t> ;
