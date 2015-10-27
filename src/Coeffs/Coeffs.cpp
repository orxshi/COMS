#include "Coeffs.h"

Coeffs::Coeffs (const Grid& gr, double rhoRef, double pRef, double MachAir, double MachAirfoil)
{
    this->rhoRef = rhoRef;
    this->pRef = pRef;
    this->MachAirfoil = MachAirfoil;
    this->MachAir = MachAir;
    totalArea = 0.;    
    dir = gr.outputDir;
    
    string tmpDir = dir;
    string temps = "lift.dat";
    string slash = "/";
    tmpDir.append (slash);
    tmpDir.append (temps);
    out.open (tmpDir);
    
    for (const Cell& cll: gr.cell)
    {
        if (cll.bc == BC::SLIP_WALL)
        {
            totalArea += mag (gr.face[cll.face[0]].area);
        }
    }
            
    double MachRef = MachAir - MachAirfoil;
    dynPresWithoutArea = 0.5 * rhoRef * pow(MachRef,2.);
    dynPresWithArea = dynPresWithoutArea * totalArea;
    
    int cntr = 0;
    for (const Cell& cll: gr.cell)
    {
        if (cll.bc == BC::SLIP_WALL)
        {
            cnt.push_back (cll.cnt);
            ++cntr;
        }
    }
    
    presCoef.resize (cntr);
}

void Coeffs::getCoeffs (const Grid& gr)
{
    CVector F;
    F[0] = 0.;
    F[1] = 0.;
    F[2] = 0.;
    
    // total force on airfoil
    int cntr = 0;
    for (const Cell& cll: gr.cell)
    {
        if (cll.bc == BC::SLIP_WALL)
        {
            presCoef[cntr] = (cll.prim[4] - pRef) / dynPresWithoutArea;
            F += cll.prim[4] * gr.face[cll.face[0]].area;
            ++cntr;
        }
    }
    
    liftCoef = F[1] / dynPresWithArea;
}

void Coeffs::outPresCoef (int id)
{
    string temps, slash, tmpDir;

    temps = "pres_";
    
    stringstream ss;
    ss << id;    
    temps.append (ss.str());
    
    temps.append (".dat");
    slash = "/";
    tmpDir = dir;
    tmpDir.append (slash);
    tmpDir.append (temps);
    ofstream outPres;
    outPres.open (tmpDir);

    if (outPres.is_open())
    {
        for (int i=0; i<cnt.size(); ++i)
        {
            outPres << setw(20) << cnt[i][0];
            outPres << setw(20) << cnt[i][1];
            outPres << setw(20) << -presCoef[i];
            outPres << endl;
        }
    }
    else
    {
        cout << "outPres is not open in outPresCoef(...)" << endl;
        exit(-2);
    }
    
    
    
}



