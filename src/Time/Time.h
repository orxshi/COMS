/* 
 * File:   Time.h
 * Author: Orhan Shibliyev
 *
 * Created on July 15, 2014, 2:37 PM
 */

#ifndef TIME_H
#define	TIME_H

#include <iostream>
#include <string>
#include <sys/time.h>

#define SEC_PER_HOUR 3600.0
#define SEC_PER_MIN 60.0

using std::cout;
using std::endl;
using std::string;

struct Watch
{
    double startTime;
    double endTime;
    double elapsedTime;
    string unit;
    bool onlySec;
    
    Watch();
    Watch (bool onlySec);
    void start();
    void stop();
    double getWallTime();    
};

#endif	/* TIME_H */

