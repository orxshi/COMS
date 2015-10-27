#include "StraightMovingAirfoil.h"

void SMAirfoil::moveEdge()
{
    double alphaRad = alpha * DEG_TO_RAD;

    CVector vel;
    vel[0] = MachAirfoil * cos (alphaRad);
    vel[1] = MachAirfoil * sin (alphaRad);
    vel[2] = 0.;

    displacePoint2D (vertex1, vel, dt, vertex2);
    displacePoint2D (vertex0, vel, dt, vertex3);
}

void SMAirfoil::read (string fileName)
{
    string tmps;
    ifstream in;
    in.open (fileName);
    
    if (in.is_open())
    {   
        in >> tmps; in >> alpha;
        in >> tmps; in >> MachAirfoil;
    }
    else
    {
        cout << "could not open file in OscAirfoil::read(...)" << endl;
        exit(-2);
    }
    
    in.close();
}

void SMAirfoil::log (string fileName)
{
    ofstream out;
    out.open (fileName, ofstream::app);
    
    if (out.is_open())
    {
        out << endl;
        out << "STRAIGHT MOVING AIRFOIL" << endl;
        
        out << "MachAirfoil = " << MachAirfoil << endl;
        out << "dt = " << dt << endl;
    }
    else
    {
        cout << "could not open file in OscAirfoil::log(...)" << endl; // log as well
        exit(-2);
    }
    
    out.close();
}