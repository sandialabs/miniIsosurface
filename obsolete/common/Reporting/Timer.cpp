/*
 * Timer.cpp
 *
 *  Created on: Jul 17, 2015
 *      Author: sjmunn
 */

#include "Timer.h"

Timer::Timer() {
	start();
}

void Timer::getCurrentTimeValues(TimeValue_t &ticksCPU,
																 TimeValue_t &wallSeconds,
																 TimeValue_t &wallNanoseconds) const {
	ticksCPU = clock();

#if (_POSIX_TIMERS > 0) && (_POSIX_MONOTONIC_CLOCK > 0)
	struct timespec currentTime;
	clock_gettime(CLOCK_MONOTONIC, &currentTime);
	wallSeconds = currentTime.tv_sec;
	wallNanoseconds = currentTime.tv_nsec;
#else
	struct timeval currentTime;
	gettimeofday(&currentTime, NULL);
	wallSeconds = currentTime.tv_sec;
	wallNanoseconds = currentTime.tv_usec * 1000;
#endif
}

void Timer::start() {
	totalWallTime=0;
	totalCPUtime=0;
	totalTicksCPU=0;

	resume();
}

void Timer::pause(void) {
	TimeValue_t currentTicksCPU, currentWallSeconds, currentWallNanoseconds;
	getCurrentTimeValues(currentTicksCPU,
											 currentWallSeconds,
											 currentWallNanoseconds);

	totalTicksCPU = currentTicksCPU-startTicksCPU;
	totalCPUtime += (static_cast<double>(totalTicksCPU))/CLOCKS_PER_SEC;

	totalWallTime += (currentWallSeconds - startWallSeconds);
	totalWallTime +=
			(currentWallNanoseconds - currentWallNanoseconds) / 1000000000.0;
}

void Timer::resume(void) {
	getCurrentTimeValues(startTicksCPU, startWallSeconds, startWallNanoseconds);
}

void Timer::stop(void) {
	pause();
}

void Timer::reportTime(YAML_Doc &doc) const {
	CLOG(logYAML) << "Total Program CPU Time (clicks): " << totalTicksCPU;
	CLOG(logYAML) << "Total Program CPU Time (seconds): " << totalCPUtime;
	CLOG(logYAML) << "Total Program WALL Time (seconds): " << totalWallTime;
	doc.add("Total Program CPU Time (clicks)",static_cast<float>(totalTicksCPU));
	doc.add("Total Program CPU Time (seconds)",totalCPUtime);
	doc.add("Total Program WALL Time (seconds)",totalWallTime);
}

Timer::~Timer() {
	// TODO Auto-generated destructor stub
}

