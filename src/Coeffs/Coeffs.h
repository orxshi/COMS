#ifndef COEFFS_H
#define	COEFFS_H

#include "../Grid/Grid.h"

struct Coeffs
{
    double rhoRef;
    double pRef;
    double MachAirfoil;
    double MachAir;
    double totalArea;
    double dynPresWithArea;
    double dynPresWithoutArea;
    double liftCoef;
    vector<double> presCoef;
    vector<CVector> cnt;
    ofstream out; // for lift    
    string dir;
    
    Coeffs (const Grid& gr, double rhoRef, double pRef, double MachAir, double MachAirfoil);
    void getCoeffs (const Grid& gr);
    void outPresCoef (int id);
};

#endif	/* COEFFS_H */

