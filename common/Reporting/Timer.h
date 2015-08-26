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

	void pause(void);
	void resume(void);

	void stop(void);
	void reportTime(YAML_Doc &) const;
	void reportTime(void) const;

private:
	struct timespec start, finish;
	clock_t tCPU;
	double totalWallTime;
	double totalCPUtime;
};

#endif /* REPORTING_TIMER_H_ */
