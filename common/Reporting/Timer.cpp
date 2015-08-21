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

void Timer::stop(void) {
	tCPU = clock()-tCPU;
	totalCPUtime=(static_cast<double>(tCPU))/CLOCKS_PER_SEC;

	clock_gettime(CLOCK_MONOTONIC, &finish);

	totalWallTime = (finish.tv_sec - start.tv_sec);
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

void Timer::reportTime(void) const {
	CLOG(logYAML) << "Total Subroutine CPU Time (clicks): " << tCPU;
	CLOG(logYAML) << "Total Subroutine CPU Time (seconds): " << totalCPUtime;
	CLOG(logYAML) << "Total Subroutine WALL Time (seconds): " << totalWallTime;
}

Timer::~Timer() {
	// TODO Auto-generated destructor stub
}

