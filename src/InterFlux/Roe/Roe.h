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

//Matrix5 jacob(const Vector<N_VAR>& q, const Vector<3>& n, double vbn);
void roeflx (Grid& gr, Limiter& limiter);

#endif	/* ROE_H */

