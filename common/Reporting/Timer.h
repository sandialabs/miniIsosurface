/*
 * Timer.h
 *
 *  Created on: Jul 17, 2015
 *      Author: sjmunn
 */

#ifndef REPORTING_TIMER_H_
#define REPORTING_TIMER_H_

#include"../includes.h"

// Reporting Headers
#include"YAML_Element.hpp"
#include"YAML_Doc.hpp"
#include"Log.h"

class Timer {
public:
	Timer();
	//Timer(bool);
	virtual ~Timer();

	//void startForProcess(void);

	void start(void);

	void pause(void);
	void resume(void);

	void stop(void);
	void reportTime(YAML_Doc &) const;
	void reportTime(void) const;

	// For message passing
	const double getWallTime (void) { return totalWallTime; };
	const double getCPUtime (void) { return totalCPUtime; };

private:
	typedef long long TimeValue_t;

	TimeValue_t startTicksCPU;
	TimeValue_t startWallSeconds;
	TimeValue_t startWallNanoseconds;
	double totalWallTime;
	double totalCPUtime;
	TimeValue_t totalTicksCPU;

	void getCurrentTimeValues(TimeValue_t &ticksCPU,
														TimeValue_t &wallSeconds,
														TimeValue_t &wallNanoseconds) const;
};

#endif /* REPORTING_TIMER_H_ */
