/*
 * MPIclockFunctions.h
 *
 *  Created on: Aug 27, 2015
 *      Author: sjmunn
 */

#ifndef MPICLOCKFUNCTIONS_H_
#define MPICLOCKFUNCTIONS_H_

// Reporting
#include"../common/Reporting/YAML_Element.hpp"
#include"../common/Reporting/YAML_Doc.hpp"
#include"../common/Reporting/Log.h"
#include"../common/Reporting/Timer.h"

#include"mpi.h"

void collectTimesMPI(YAML_Doc &doc, Timer * runTime) {
	//  Get the individual process ID.
	int pID = MPI::COMM_WORLD.Get_rank();
	int nProcesses = MPI::COMM_WORLD.Get_size();

	double wallTime = runTime->getWallTime();
	double *wallTimes = NULL;
	if (pID == 0) {
		wallTimes = new double [nProcesses];
	}
	MPI_Gather(&wallTime, 1, MPI_DOUBLE, wallTimes, 1, MPI_DOUBLE, 0,
	           MPI_COMM_WORLD);
	double CPUtime = runTime->getWallTime();
	double *CPUtimes = NULL;
	if (pID == 0) {
		CPUtimes = new double [nProcesses];
	}
	MPI_Gather(&CPUtime, 1, MPI_DOUBLE, CPUtimes, 1, MPI_DOUBLE, 0,
	           MPI_COMM_WORLD);

	// Generate a final report of the execution times
	if (pID == 0) {
		double maxWallTime=0;
		double totalCPU=0;
		for (int iPID=0;iPID<nProcesses;iPID++) {
			CLOG(logYAML) << "Process number, " << iPID << " CPU Time (seconds): " << CPUtimes[iPID];
			CLOG(logYAML) << "Process number, " << iPID << " Wall Time (seconds): " << wallTimes[iPID];
			std::string reportCPU = "Process number, "
					+ std::to_string(static_cast<long long int>(iPID)) + " CPU Time (seconds)";
			std::string reportWall = "Process number, "
					+ std::to_string(static_cast<long long int>(iPID)) + " Wall Time (seconds):";
			doc.add(reportCPU,wallTimes[iPID]);
			doc.add(reportWall,wallTimes[iPID]);
			if (wallTimes[iPID] > maxWallTime) maxWallTime=wallTimes[iPID];
			totalCPU+=CPUtimes[iPID];
		}
		CLOG(logYAML) << "Total CPU time: " << totalCPU;
		CLOG(logYAML) << "Overall wall time: " << maxWallTime;
		doc.add("Total CPU time: ",totalCPU);
		doc.add("Overall wall time: ",maxWallTime);
	}
	doc.generateYAML();
}

#endif /* MPICLOCKFUNCTIONS_H_ */
