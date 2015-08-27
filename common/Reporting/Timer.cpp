/*
 * Timer.cpp
 *
 *  Created on: Jul 17, 2015
 *      Author: sjmunn
 */

#include "Timer.h"

Timer::Timer() {
	totalWallTime=0;
	totalCPUtime=0;
	tCPU = clock();
	clock_gettime(CLOCK_MONOTONIC, &start);
}

//Timer::Timer(bool hold) {
//	totalWallTime=0;
//	totalCPUtime=0;
//	tCPU = 0;
//}
//
//void Timer::startForProcess(void) {
//	totalWallTime=0;
//	totalCPUtime=0;
//	tCPU = clock();
//	clock_gettime(CLOCK_MONOTONIC, &start);
//}

void Timer::pause(void) {
	tCPU = clock()-tCPU;
	totalCPUtime += (static_cast<double>(tCPU))/CLOCKS_PER_SEC;

	clock_gettime(CLOCK_MONOTONIC, &finish);

	totalWallTime += (finish.tv_sec - start.tv_sec);
	totalWallTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
}

void Timer::resume(void) {
	tCPU = clock();
	clock_gettime(CLOCK_MONOTONIC, &start);
}

void Timer::stop(void) {
	tCPU = clock()-tCPU;
	totalCPUtime += (static_cast<double>(tCPU))/CLOCKS_PER_SEC;

	clock_gettime(CLOCK_MONOTONIC, &finish);

	totalWallTime += (finish.tv_sec - start.tv_sec);
	totalWallTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
}

void Timer::reportTime(YAML_Doc &doc) const {
	CLOG(logYAML) << "Total Program CPU Time (clicks): " << tCPU;
	CLOG(logYAML) << "Total Program CPU Time (seconds): " << totalCPUtime;
	CLOG(logYAML) << "Total Program WALL Time (seconds): " << totalWallTime;
	doc.add("Total Program CPU Time (clicks)",static_cast<float>(tCPU));
	doc.add("Total Program CPU Time (seconds)",totalCPUtime);
	doc.add("Total Program WALL Time (seconds)",totalWallTime);
}

//void Timer::reportTimeMPI(YAML_Doc &doc, int pID) const {
//	// Reporting process times to the master thread
//	int position = 0;
//	char buffer[20];
//	if (pID ==1) {
//		// Senders
//		const double wallTime = runTime.getWallTime();
//		const double CPUtime = runTime.getCPUtime();
//
//		MPI_Pack(&wallTime, 1, MPI_DOUBLE, buffer, 20, &position, MPI_COMM_WORLD);
//		MPI_Pack(&CPUtime, 1, MPI_DOUBLE, buffer, 20, &position, MPI_COMM_WORLD);
//
//		MPI_Send(buffer, position, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
//	}
//	else if (pID == 0) {
//		double wallTime=0,CPUtime=0;
//
//		MPI_Status status;
//		// Receiver
//		CLOG(logDEBUG) << "Waiting for message";
//		MPI_Recv(buffer, 110, MPI_PACKED, 1, 0, MPI_COMM_WORLD, &status);
//		CLOG(logDEBUG) << "Message received";
//		MPI_Unpack(buffer,1000,&position,&wallTime,1,MPI_DOUBLE,MPI_COMM_WORLD);
//		MPI_Unpack(buffer,1000,&position,&CPUtime,1,MPI_DOUBLE,MPI_COMM_WORLD);
//
//		CLOG(logDEBUG) << "Received Program CPU Time (seconds): " << CPUtime;
//		CLOG(logDEBUG) << "Received Program WALL Time (seconds): " << wallTime;
//	}
//
//	CLOG(logYAML) << "Total Program CPU Time (clicks): " << tCPU;
//	CLOG(logYAML) << "Total Program CPU Time (seconds): " << totalCPUtime;
//	CLOG(logYAML) << "Total Program WALL Time (seconds): " << totalWallTime;
//	doc.add("Total Program CPU Time (clicks)",static_cast<float>(tCPU));
//	doc.add("Total Program CPU Time (seconds)",totalCPUtime);
//	doc.add("Total Program WALL Time (seconds)",totalWallTime);
//}

Timer::~Timer() {
	// TODO Auto-generated destructor stub
}

