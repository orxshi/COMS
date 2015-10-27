#include "Solid.h"

void Solid::setSurfaceArea()
{
    surfaceArea = 0.;
    for (const CVector& a: area)
    {
        surfaceArea += mag(a);
    }
}

void Solid::setDynPressure()
{
    dynP = 0.5 * rhoRef * pow(machRef,2) * surfaceArea;
}

void Solid::getPC()
{
    if (pc.empty())
    {
        pc.resize(pres.size());
    }
    
    double tmp = ( GAMMA * pow(machRef, 2) );
    
    for (int i=0; i<pres.size(); ++i)
    {
        pc[i] = -2. * (pres[i] / pRef - 1.) / tmp;
    }
}

void Solid::outPC (string dir)
{
    string temps = "pc.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);

    ofstream out;
    out.open (dir);
    
    for (int i=0; i<pc.size(); ++i)
    {
        out << setw(20) << cnt[i][0];
        out << setw(20) << pc[i];
        out << endl;
    }
    
    out.close();
}

void Solid::outLDM (string dir)
{
    string temps = "ldm.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);

    ofstream out;
    out.open (dir);
    
    out << setw(20) << dc;
    out << setw(20) << lc;
    out << setw(20) << mc;
    out << endl;
    
    out.close();
}

void Solid::getLDM()
{
    CVector force;
    CVector moment;
    CVector r, n1, fn;    
    //string dir = outputDir;
    
    //double relMach = MachAir - MachAirfoil;

    force.fill(0.);
    moment.fill(0.);

    // force & moment on airfoil surface        
    for (int i=0; i<pres.size(); ++i)
    {
        force += ( ( pres[i] - pRef ) * area[i] );

        r = cnt[i] - solidCnt;
        r[2] = 0.;
        n1[0] = -r[1];
        n1[1] =  r[0];
        n1[2] = 0.;
        fn = dotP( force,norm(n1) ) * norm(n1);
        moment += crossP( fn,r );
    }

    double Fd = force[0];
    double Fl = force[1];

    dc = Fd / dynP; // drag coef
    lc = Fl / dynP; // lift coef
    mc = mag(moment) / ( dynP / surfaceArea ); // moment coef
}

void Solid::read()
{
    ifstream in;
    in.open ("solidData.dat");
    
    in >> phys;
    in >> chord;
    in >> solidCnt[0];
    in >> solidCnt[1];
    in >> solidCnt[2];
    
    in.close();
}