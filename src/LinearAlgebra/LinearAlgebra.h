/* 
 * File:   LinearAlgebra.h
 * Author: Orhan Shibliyev
 *
 * Created on July 22, 2014, 11:49 AM
 */

#ifndef LINEARALGEBRA_H
#define	LINEARALGEBRA_H

#include <vector>
#include "../Vector/Vector.h"
#include "../Constants.h"

using std::vector;

void GaussElim (Vector2D<N_DIM,N_DIM> A, CVector b, CVector& x);
void newtonSolve(const Vector2D<N_DIM,8>& f, double& u, double& v, double& w);
void osInterpolants (const vector<CVector>& xc, const CVector& xp, const unsigned int nVtx, double frac[]);

#endif	/* LINEARALGEBRA_H */

