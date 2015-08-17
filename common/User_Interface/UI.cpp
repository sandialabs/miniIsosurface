/*
 * UI.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: sjmunn
 */

#include "UI.h"

UI::UI(int argc, char* argv[]) {

	filePath = argv[1];
	outFilePath = argv[2];
	isoval = atof(argv[3]);
	this->checkDebugLevel();
}

const char * UI::getFile(void) const {
	return filePath;
}

const char * UI::outFile(void) const {
	return outFilePath;
}

float_t UI::getIsoVal(void) const {
	return isoval;
}

void UI::checkDebugLevel(void) const {
	if (LOG::ReportingLevel() == logDEBUG_Step) {
		CLOG(logWARNING) << "Step by step debugging output is active";
		CLOG(logWARNING)
				<< "Only image volumes with up to 1000 points will be accepted";
	}
}

void UI::checkArgs(int argc) {
	if (argc != 4) {
		throw wrong_arguments(argc);
	}
}

//void UI::marchImplemtation(const Image3D_t & vol, TriangleMesh_t *&mesh,
//		YAML_Doc &doc) {
//	mergemp::extractIsosurface(vol, isoval, mesh, doc);
//}
