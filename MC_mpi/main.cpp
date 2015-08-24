/*
 * main.cpp
 *
 *  Created on: Jul 2, 2015
 *      Author: sjmunn
 */

// Common utility headers -------------
// Standard C/C++ library
#include"../common/includes.h"
#include"mpi.h"

// File i/o
#include"../common/IO/LoadImage3DMPI.h"
#include"../common/IO/SaveTriangleMesh.h"

// Reporting
#include"../common/Reporting/YAML_Element.hpp"
#include"../common/Reporting/YAML_Doc.hpp"
#include"../common/Reporting/Log.h"
#include"../common/Reporting/Timer.h"

#include"../common/User_Interface/UI.h"

// Data Objects
#include"../common/types.h"

// Local implementation headers -------
// Algorithm
#include"./Implementations/MpiAlgo.h"

int main(int argc, char* argv[]) {

	LOG::ReportingLevel() = logDEBUG1; // Debug level is hard coded

	// Initialize the user interface
	UI<float_t> mainUI(argc,argv);

	// Create data object
	GeneralContext<float_t> data;

	// Only load the header file, MpiAlgo will read the rest
	LoadImage3DMPI<float_t> fileHeader;
	fileHeader.loadHeader(mainUI.getFile());

	// Report file data characteristics
	CLOG(logYAML) << "Volume image data file path: " << mainUI.getFile();
	data.doc.add("Volume image data file path", mainUI.getFile());
	fileHeader.report(data.doc);

	// Initialize MPI
	MPI::Init(argc, argv);
	// Start the clock
	Timer RunTime;
	// Execute the marching cubes implementation
	data.isoval=mainUI.getIsoVal();
	MpiAlgo<float_t> algorithm(fileHeader); // Specific to MPI
	data.setAlgorithm(&algorithm);
	data.march();
	// Stop Clock
	RunTime.stop();

	//Report and save YAML file
	RunTime.reportTime(data.doc);
	data.doc.generateYAML();

	// Save the result
	saveTriangleMesh(&(data.mesh), mainUI.outFile());
	return 0;
}
