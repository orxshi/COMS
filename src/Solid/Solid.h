/* 
 * File:   Solid.h
 * Author: Orhan Shibliyev
 *
 * Created on October 6, 2014, 12:02 AM
 */

#ifndef SOLID_H
#define	SOLID_H

#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../Constants.h"
#include "../Vector/Vector.h"

using std::pow;
using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::setw;
using std::endl;

struct Solid
{
    int phys;
    double dynP;
    double rhoRef;
    double pRef;
    double machRef;
    double surfaceArea;
    double chord;
    double lc, dc, mc;    
    CVector solidCnt;
    vector<double> pc;
    vector <double> pres;
    vector <CVector> area;
    vector <CVector> cnt;
  
    void setSurfaceArea();
    void setDynPressure();
    void getPC();
    void getLDM();
    void outPC (string dir);
    void outLDM (string dir);
    void read();
};

#endif	/* SOLID_H */

