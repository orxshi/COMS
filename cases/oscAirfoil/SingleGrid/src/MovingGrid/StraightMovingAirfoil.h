/* 
 * File:   StraightMovingAirfoil.h
 * Author: orhan
 *
 * Created on June 27, 2015, 12:34 PM
 */

#ifndef STRAIGHTMOVINGAIRFOIL_H
#define	STRAIGHTMOVINGAIRFOIL_H

#include "../../../../../src/MovingMesh/MovingMesh.h"

struct SMAirfoil:MovingGrid
{
    double alpha;    
    double MachAirfoil;
    
    SMAirfoil (double dt) : MovingGrid (dt) {}
    
    void moveEdge();
    void moveGrid (Grid& gr) {} // no need to move
    void read (string fileName);
    void log (string fileName);
    void interFromOldTS (Grid& curGrid, Grid& oldGrid) {} // not defined yet
};

#endif	/* STRAIGHTMOVINGAIRFOIL_H */

