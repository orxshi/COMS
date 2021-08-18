#ifndef OSCINIT_H
#define	OSCINIT_H

#include "../../../../../src/Init/Init.h"

struct OscInit:Init
{
    double rhoInf;
    double Mach;
    double aoa;
    double pInf;
   	
    void read();
    void init (Grid& gr);
    void init_sod (Grid& gr);
    void print();
    void log (string dir);
};

#endif
