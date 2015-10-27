#include "Time.h"

Watch::Watch()
{
    unit = "sec";
    onlySec = false;
}

Watch::Watch (bool onlySec)
{
    this->onlySec = onlySec;
}

double Watch::getWallTime()
{
    struct timeval time;
    
    gettimeofday(&time,NULL);
    
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void Watch::start()
{
    startTime = getWallTime();
}

void Watch::stop()
{
    endTime = getWallTime();
    
    elapsedTime = endTime - startTime;
    
    if (!onlySec)
    {
        if (elapsedTime > SEC_PER_HOUR)
        {
            elapsedTime /= SEC_PER_HOUR;
            unit = "hour";
        }
        else if (elapsedTime > SEC_PER_MIN)
        {
            elapsedTime /= SEC_PER_MIN;
            unit = "min";
        }
    }
}