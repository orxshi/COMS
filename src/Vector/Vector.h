/* 
 * File:   Vector.h
 * Author: Orhan Shibliyev
 *
 * Created on June 17, 2014, 1:30 AM
 */

#ifndef VECTOR_H
#define	VECTOR_H

#include <array>
#include <cmath>
#include "../Constants.h"

using std::array;
using std::size_t;
using std::sqrt;

template<size_t s> using Vector = array <double, s>;
template<size_t sI, size_t sO> using Vector2D = array <Vector<sI>, sO>;
using CVector = Vector<N_DIM>; // N-dimensional Cartesian vector

// Operator Overloading
template<size_t s> Vector<s> operator+ (Vector<s> V1, const Vector<s>& V2)
{
    for (int i=0; i<s; ++i)
    {
        V1[i] += V2[i];
    }
    
    return V1;
}
template<size_t s> Vector<s> operator- (Vector<s> V1, const Vector<s>& V2)
{
    for (int i=0; i<s; ++i)
    {
        V1[i] -= V2[i];
    }
    
    return V1;
}
template<size_t s> Vector<s> operator* (const double sc, Vector<s> V)
{
    for (int i=0; i<s; ++i)
    {
        V[i] *= sc;
    }

    return V;
}
template<size_t s> void operator+= (Vector<s>& V1, const Vector<s>& V2)
{
    for (int i=0; i<s; ++i)
    {
        V1[i] += V2[i];
    }
}
template<size_t s> void operator/= (Vector<s>& V, const double scalar)
{
    for (int i=0; i<s; ++i)
    {
        V[i] /= scalar;
    }
}
template<size_t s> void operator-= (Vector<s>& V1, const Vector<s>& V2)
{
    for (int i=0; i<s; ++i)
    {
        V1[i] -= V2[i];
    }
}
template<size_t s> void operator*= (Vector<s>& V, const double scalar)
{
    for (int i=0; i<s; ++i)
    {
        V[i] *= scalar;
    }
}
template<size_t s> double mag(const Vector<s>& V)
{
    double mg = 0.;
    
    for (int i=0; i<s; ++i)
    {
        mg += pow(V[i],2);
    }
    
    mg = sqrt(mg);
    
    return mg;
}
template<size_t s> Vector<s> norm(Vector<s> V)
{
    double mg = mag(V);
    
    for (int i=0; i<s; ++i)
    {
        V[i] /= mg;
    }
    
    return V;
}
double dotP (const Vector<3>& V1, const Vector<3>& V2);
Vector<3> crossP (const Vector<3>& V1, const Vector<3>& V2);

#endif	/* VECTOR_H */

