/* 
 * File:   Roe.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 7:57 AM
 */

#ifndef ROE_H
#define	ROE_H

#include <array>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include "../../Constants.h"
#include "../../Vector/Vector.h"
#include "../../Matrix5/Matrix5.h"
#include "../../Face/Face.h"
#include "../../Grid/Grid.h"
#include "../../Limiter/Limiter.h"


using std::array;
using std::pow;
using std::cout;
using std::endl;

struct Roe
{
    int rank;
    int nProcs;
    
    int localSizeFace;
    int* localSizesFace;
    int* displsFace;
    
    // M
    int localSizeFaceM;
    int* localSizesFaceM;
    int* displsFaceM;
    
    // vel
    int localSizeFaceV;
    int* localSizesFaceV;
    int* displsFaceV;
    
    // flux
    int localSizeFaceF;
    int* localSizesFaceF;
    int* displsFaceF;
    
    vector <double> vel;
    vector <Vector<N_VAR>> flux;
    
    Roe (Grid& gr);
    void initParallelVars (Grid& gr);
    void roeflx (Grid& gr, Limiter& limiter, vector <Matrixd<N_VAR,N_VAR>>& M0, vector <Matrixd<N_VAR,N_VAR>>& M1);
};

//Matrix5 jacob(const Vector<N_VAR>& q, const Vector<3>& n, double vbn);


#endif	/* ROE_H */

