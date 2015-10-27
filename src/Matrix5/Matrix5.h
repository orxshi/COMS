/* 
 * File:   Matrix5.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 8:00 AM
 */

#ifndef MATRIX5_H
#define	MATRIX5_H

#include <array>
#include <iostream>
#include "../Vector/Vector.h"

using std::array;
using std::cout;
using std::endl;

#define MAT5_SIZE 5

using Matrix5 = array < array<double, MAT5_SIZE>, MAT5_SIZE >;

Matrix5 add5 (const Matrix5& M1, const Matrix5& M2);
Matrix5 sub5 (const Matrix5& M1, const Matrix5& M2);
Matrix5 mul5 (const Matrix5& M1, const Matrix5& M2);
Matrix5 mulI5 (const Matrix5& M1, const array<double,MAT5_SIZE>& a);
Matrix5 mulS5 (const double s, const Matrix5& M);
void eq5 (Matrix5& M, const double s);
Vector<MAT5_SIZE> mat5Vec5Mul (const Matrix5& M, const Vector<MAT5_SIZE>& V);

#endif	/* MATRIX5_H */
