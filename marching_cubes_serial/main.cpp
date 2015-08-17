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
#include"../common/IO/LoadImage3D.h"
#include"../common/IO/SaveTriangleMesh.h"

// Reporting
#include"../common/Reporting/YAML_Element.hpp"
#include"../common/Reporting/YAML_Doc.hpp"
#include"../common/Reporting/Log.h"
#include"../common/Reporting/Timer.h"

// Data Objects
#include"../common/types.h"

// Local implementation headers -------
// User interface
#include"./User_Interface/SerialInterface.h"

// Implementation Data
#include"./Implementations/SerialAlgo.h"

typedef SerialData<float_t> SerialData_t;

int main(int argc, char* argv[]) {

	LOG::ReportingLevel() = logDEBUG; // Debug level is hard coded

	// Initialize the user interface
	SerialInterface<float_t> mainUI(argc,argv);

	// Initialize console log and YAML
	YAML_Doc doc("Marching Cubes", "0.1", ".", "yaml_out.yaml");

	// Create runtime object
	SerialData_t runData;
	loadImage3D(mainUI.getFile(), &(runData.imageIn));

	// Report file data characteristics
	CLOG(logYAML) << "Volume image data file path: " << mainUI.getFile();
	doc.add("Volume image data file path", mainUI.getFile());
	runData.imageIn.report(doc);

	// Start the clock
	Timer RunTime;

	mainUI.marchImplemtation(runData);

	// Stop Clock
	RunTime.stop();
	//Report and save YAML file
	RunTime.reportTime(doc);
	doc.generateYAML();

	// Save the result
	saveTriangleMesh(&(runData.mesh), mainUI.outFile());
	return 0;
}
