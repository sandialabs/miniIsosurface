/*
 * main.cpp
 *
 *  Created on: Jul 2, 2015
 *      Author: sjmunn
 */

// Common utility headers -------------
// Standard C/C++ library
#include"../common/includes.h"

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

// Deals with each thread's timer:
#include"./MPIclockFunctions.h"

int main(int argc, char* argv[]) {

	LOG::ReportingLevel() = logDEBUG1; // Debug level is hard coded

	int id;
	int nProcesses;

	// Initialize MPI
	MPI::Init(argc, argv);
	//  Get the number of processes.
	nProcesses = MPI::COMM_WORLD.Get_size();
	//  Get the individual process ID.
	id = MPI::COMM_WORLD.Get_rank();

	// Initialize the user interface
	UI<float_t> mainUI(argc, argv);

	// Create data object
	GeneralContext<float_t> data;

	// Only load the header file, MpiAlgo will read the rest
	LoadImage3DMPI<float_t> fileHeader;
	fileHeader.loadHeader(mainUI.getFile());

	// Report file data characteristics
	CLOG(logYAML) << "Volume image data file path: " << mainUI.getFile();
	data.doc.add("Volume image data file path", mainUI.getFile());
	fileHeader.report(data.doc);

	// Start the clock
	Timer runTime;
	// Execute the marching cubes implementation
	data.isoval = mainUI.getIsoVal();
	MpiAlgo<float_t> algorithm(fileHeader, id, nProcesses, &runTime);
	data.setAlgorithm(&algorithm);
	data.march();
	// Stop Clock
	runTime.stop();

	collectTimesMPI(data.doc,&runTime);

	// Save the result
	const char * baseName = mainUI.outFile();
	std::string outFile = baseName;
	outFile = outFile + "." + std::to_string(static_cast<long long int>(id));
	saveTriangleMesh(&(data.mesh), outFile.c_str());

	//
	//  Terminate MPI.
	//
	MPI::Finalize();
	return 0;
}
