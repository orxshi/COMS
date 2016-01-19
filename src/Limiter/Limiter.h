/* 
 * File:   Limiter.h
 * Author: Orhan Shibliyev
 *
 * Created on July 22, 2014, 6:41 PM
 */

#ifndef LIMITER_H
#define	LIMITER_H

#include "../Vector/Vector.h"
#include "../Constants.h"
//#include "../Grid/Grid.h"
#include "../Gradient/Gradient.h"
#include <algorithm>
#include <vector>
#include <mpi.h>

using std::min;
using std::max;
using std::pow;
using std::vector;

//void minMod(const Vector2D<3,N_VAR>& gradL, const Vector2D<3,N_VAR>& gradR, Vector2D<3,N_VAR>& grad);
//double venkata (const Vector2D<3,N_VAR>& grad, Vector<N_VAR>& uMax, Vector<N_VAR>& uMin, Vector<N_VAR>& u, const CVector& dis);

struct Limiter
{
    int type;
    // minmod = 0
    // Van Albada = 1
    // bj = 2
    // venka = 3
    vector<vector<double>> ksiBJ;
    vector<vector<double>> ksiV;
    
    int rank;
    int nProcs;
    int localSize;
    int localSizeNVAR;
    int* localSizes;    
    int* localSizesNVAR;
    int* displs;
    int* displsNVAR;
    
    Limiter (Grid& gr);
    
    void getLimitedGrad (const Vector2D<3,N_VAR>& gradL, const Vector2D<3,N_VAR>& gradR, Vector2D<3,N_VAR>& grad); // not verified
    void getLimitedGradDarwish (const Vector2D<3,N_VAR>& gradL, const Vector2D<3,N_VAR>& gradR, const CVector& rL, const CVector& rR,
                                const Vector<N_VAR>& varL, const Vector<N_VAR>& varR, Vector<N_VAR>& reconstL, Vector<N_VAR>& reconstR); // includes minmod only
    void bj (Grid& gr, Gradient& gradient);
    void venka (Grid& gr, Gradient& gradient);
    void initParallelVars (Grid& gr);
};

#endif	/* LIMITER_H */

