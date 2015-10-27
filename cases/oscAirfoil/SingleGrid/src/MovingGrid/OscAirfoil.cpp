#include "OscAirfoil.h"

void OscAirfoil::setAngles(const double time)
{
    double tmp = 2. * kc * fabs(MachAirfoil);
    //alpha    = alphaMean + alphaMax * sin( tmp * (time+dt) ); // degree
    //delAlpha = -dt * tmp * alphaMax * cos( tmp * (time+dt) ); // degree
    
    alpha    = alphaMean + alphaMax * sin( tmp * time ); // degree
    delAlpha = -tmp * alphaMax * cos( tmp * time ); // degree
}

void OscAirfoil::moveEdge()
{
    CVector r;
    
    double delAlphaRad = delAlpha * DEG_TO_RAD;
    double alphaRad = alpha * DEG_TO_RAD;
    
    r = vertex1 - centerAirfoil;
    rotateVectorAroundPoint2D (centerAirfoil, delAlphaRad, r, vertex2);

    r = vertex0 - centerAirfoil;
    rotateVectorAroundPoint2D (centerAirfoil, delAlphaRad, r, vertex3);

    CVector vel;
    vel[0] = MachAirfoil * cos (alphaRad);
    vel[1] = MachAirfoil * sin (alphaRad);
    vel[2] = 0.;

    displacePoint2D (vertex2, vel, dt, vertex2);
    displacePoint2D (vertex3, vel, dt, vertex3);
}

void OscAirfoil::moveGrid (Grid& gr)
{
    CVector r;

    double dar = delAlpha * DEG_TO_RAD;
    
    for (Point& p: gr.pt)
    {
        r = p.dim - centerAirfoil;
        rotateVectorAroundPoint2D (centerAirfoil, dar, r, p.dim);
    }
    
    for (Cell& c: gr.cell) { c.set_centroid (gr.pt); }

    for (Face& f: gr.face)
    {
        f.set_centroid (gr.pt);
        f.set_area (gr.pt);
    }

    for (Cell& c: gr.cell)
    {
        double oldU = c.prim[1];
        double oldV = c.prim[2];

        c.prim[1] = oldU * cos(dar) - oldV * sin(dar);
        c.prim[2] = oldU * sin(dar) + oldV * cos(dar);
        
        c.prim_to_cons();
    }
}

void OscAirfoil::read (string fileName)
{
    string tmps;
    ifstream in;
    in.open (fileName);
    
    if (in.is_open())
    {
        in >> tmps; in >> alphaMean;
        in >> tmps; in >> alphaMax;
        in >> tmps; in >> kc;
        in >> tmps; in >> MachAirfoil;
        in >> tmps; in >> centerAirfoil[0];
        in >> tmps; in >> centerAirfoil[1];
        in >> tmps; in >> centerAirfoil[2];
    }
    else
    {
        cout << "could not open file in OscAirfoil::read(...)" << endl;
        exit(-2);
    }
    
    in.close();
}

void OscAirfoil::log (string fileName)
{
    ofstream out;
    out.open (fileName, ofstream::app);
    
    if (out.is_open())
    {
        out << endl;
        out << "OSCILLATING AIRFOIL" << endl;
        
        out << "alphaMean = " << alphaMean << endl;
        out << "alphaMax = " << alphaMax << endl;
        out << "kc = " << kc << endl;
        out << "MachAirfoil = " << MachAirfoil << endl;
        out << "centerAirfoil[0] = " << centerAirfoil[0] << endl;
        out << "centerAirfoil[1] = " << centerAirfoil[1] << endl;
        out << "centerAirfoil[2] = " << centerAirfoil[2] << endl;
        out << "dt = " << dt << endl;
    }
    else
    {
        cout << "could not open file in OscAirfoil::log(...)" << endl; // log as well
        exit(-2);
    }
    
    out.close();
}