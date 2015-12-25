/* 
 * File:   Matrix5.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 8:00 AM
 */

#ifndef MATRIX5_H
#define	MATRIX5_H

#include <array>
#include <vector>
#include <iostream>
#include "../Vector/Vector.h"

using std::array;
using std::cout;
using std::endl;
using std::vector;

#define MAT5_SIZE 5

using Matrix5 = array < array<double, MAT5_SIZE>, MAT5_SIZE >;

Matrix5 add5 (const Matrix5& M1, const Matrix5& M2);
Matrix5 sub5 (const Matrix5& M1, const Matrix5& M2);
Matrix5 mul5 (const Matrix5& M1, const Matrix5& M2);
Matrix5 mulI5 (const Matrix5& M1, const array<double,MAT5_SIZE>& a);
Matrix5 mulS5 (const double s, const Matrix5& M);
void eq5 (Matrix5& M, const double s);
Vector<MAT5_SIZE> mat5Vec5Mul (const Matrix5& M, const Vector<MAT5_SIZE>& V);

template<size_t nRow, size_t nCol> struct Matrixd
{
    //int dim;
    //vector<int> sizes;
    array<double, nRow * nCol> mat;
    
    double operator()(const int row, const int col) const
    {
        return mat[row*nCol + col];
    }

    double& operator()(const int row, const int col)
    {
        return mat[row*nCol + col];
    }
    
    template<size_t nRowO, size_t nColO> Matrixd<nRow, nColO> operator* (const Matrixd<nRowO, nColO>& M) const
    {
        // matrix-matrix product
    
        if (nRowO != nCol)
        {
            cout << "mismatch of sizes in Matrixd::mulMM(...)" << endl;
            exit(-2);
        }
    
        Matrixd<nRow, nColO> T;
        
        for (int r=0; r<nRow; ++r)
        {
            for (int c=0; c<nColO; ++c)
            {
                T(r,c) = 0.; // take a row
                
                for (int j=0; j<nCol; ++j)
                {
                    T(r,c) += mat[r*nCol+j] * M(j,c);
                }
            }
        }
        
        return T;
    }
    
    Matrixd operator* (const array<double,nCol>& v) const
    {
        // matrix-diag matrix product
        // v is nCol x nCol diagonal square matrix
        
        Matrixd<nRow, nCol> T;

        for (int r=0; r<nRow; ++r)
        {
            for (int c=0; c<nCol; ++c)
            {
                T(r,c) = mat[r*nCol+c] * v[c];
            }
        }
        
        return T;
    }
    
    Vector<nRow> operator% (const array<double,nCol>& v) const
    {
        // matrix-vector product
    
        Vector<nRow> V;
        V.fill(0.);
    
        for (int r=0; r<nRow; ++r)
        {
            for (int c=0; c<nCol; ++c)
            {
                V[r] += mat[r*nCol+c] * v[c];
            }
        }
        
        return V;
    }
    
    Matrixd operator* (const double s) const
    {
        Matrixd<nRow, nCol> T;

        for (int r=0; r<nRow; ++r)
        {
            for (int c=0; c<nCol; ++c)
            {
                T(r,c) = s * mat[r*nCol+c];
            }
        }
        
        return T;
    }
    
    Matrixd operator+ (const Matrixd<nRow, nCol>& M) const
    {
        Matrixd<nRow, nCol> T;

        for (int r=0; r<nRow; ++r)
        {
            for (int c=0; c<nCol; ++c)
            {
                T(r,c) = mat[r*nCol+c] + M(r,c);
            }
        }
        
        return T;
    }
    
    Matrixd operator- (const Matrixd<nRow, nCol>& M) const
    {
        Matrixd<nRow, nCol> T;

        for (int r=0; r<nRow; ++r)
        {
            for (int c=0; c<nCol; ++c)
            {
                T(r,c) = mat[r*nCol+c] - M(r,c);
            }
        }
        
        return T;
    }
    
    Matrixd& operator= (const double s)
    {
        for (int r=0; r<nRow; ++r)
        {
            for (int c=0; c<nCol; ++c)
            {
                mat[r*nCol+c] = s;
            }
        }
    }
};

#endif	/* MATRIX5_H */
