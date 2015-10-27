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

void minMod(const Vector2D<3,N_VAR>& gradL, const Vector2D<3,N_VAR>& gradR, Vector2D<3,N_VAR>& grad);

#endif	/* LIMITER_H */

