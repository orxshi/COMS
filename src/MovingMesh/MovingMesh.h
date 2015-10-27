/* 
 * File:   MovingMesh.h
 * Author: Orhan Shibliyev
 *
 * Created on September 22, 2014, 5:59 PM
 */

#ifndef MOVINGMESH_H
#define	MOVINGMESH_H

#include "../Grid/Grid.h"

struct MovingGrid
{
    CVector vertex0, vertex1, vertex2, vertex3;
    double dt;
    
    MovingGrid (double dt);
    
    int getEdge (const vector<Point>& pt, const Face& f);
    virtual void moveEdge() = 0;
    void getFaceVelocity(Face& f);
    void checkGCL(const Grid& gr);
    void getAllFaceVelocities (Grid& gr);
    virtual void moveGrid (Grid& gr) = 0;
    virtual void read (string fileName) = 0;
    void rotateVectorAroundPoint2D (const CVector pivot, const double angleRad, const CVector r, CVector& ans);
    void displacePoint2D (const CVector orig, const CVector vel, const double dt, CVector& ans);
    virtual void log (string fileName) = 0;
    virtual void interFromOldTS (Grid& curGrid, Grid& oldGrid) = 0;
};

#endif	/* MOVINGMESH_H */

