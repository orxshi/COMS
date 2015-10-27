/* 
 * File:   OscAirfoil.h
 * Author: orhan
 *
 * Created on June 2, 2015, 8:24 PM
 */

#ifndef OSCAIRFOIL_H
#define	OSCAIRFOIL_H

#include "../../../../../src/MovingMesh/MovingMesh.h"

struct OscAirfoil:MovingGrid
{
    double alpha;
    double delAlpha;
    double alphaMean;
    double alphaMax;
    double kc;
    double MachAirfoil;    
    //double dt;
    //double time;
    CVector centerAirfoil;
    
    OscAirfoil (double dt) : MovingGrid (dt) {}
    
    void moveEdge();
    void moveGrid (Grid& gr);
    void setAngles (const double time);
    void read (string fileName);
    void log (string fileName);
    void interFromOldTS (Grid& curGrid, Grid& oldGrid);
};

#endif	/* OSCAIRFOIL_H */

