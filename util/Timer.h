/*
 * Timer.h
 *
 *  Created on: Jul 17, 2015
 *      Author: sjmunn
 */

#ifndef REPORTING_TIMER_H_
#define REPORTING_TIMER_H_

#include <stdlib.h>
#include <ctime>
#include <sys/time.h>

namespace util {

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

    double getTotalTicks(void) const { return totalTicksCPU; }
    double getWallTime(void) const   { return totalWallTime; }
    double getCPUtime(void) const    { return totalCPUtime; }

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

}

#endif /* REPORTING_TIMER_H_ */
