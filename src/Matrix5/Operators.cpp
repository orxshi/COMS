#include "Matrix5.h"

Matrix5 add5 (const Matrix5& M1, const Matrix5& M2)
{
    Matrix5 T;

    for (int r=0; r<MAT5_SIZE; ++r)
    {
        for (int c=0; c<MAT5_SIZE; ++c)
        {
            T[r][c] = M1[r][c] + M2[r][c];
        }
    }
    
    return T;
}

Matrix5 sub5 (const Matrix5& M1, const Matrix5& M2)
{
    Matrix5 T;

    for (int r=0; r<MAT5_SIZE; ++r)
    {
        for (int c=0; c<MAT5_SIZE; ++c)
        {
            T[r][c] = M1[r][c] - M2[r][c];
        }
    }
    
    return T;
}

Matrix5 mulI5 (const Matrix5& M1, const array<double,MAT5_SIZE>& a)
{
    Matrix5 T;

    for (int r=0; r<MAT5_SIZE; ++r)
    {
        for (int c=0; c<MAT5_SIZE; ++c)
        {
            T[r][c] = M1[r][c] * a[c];
        }
    }
    
    return T;
}

Matrix5 mul5 (const Matrix5& M1, const Matrix5& M2)
{
    Matrix5 T;
    
    for (int r=0; r<MAT5_SIZE; ++r)
    {
        for (int c=0; c<MAT5_SIZE; ++c)
        {
            T[r][c] = 0.; // take a row
            
            for (int j=0; j<MAT5_SIZE; ++j)
            {
                T[r][c] += M1[r][j] * M2[j][c];
            }
        }
    }
    
    return T;
}

Matrix5 mulS5 (const double s, const Matrix5& M)
{
    Matrix5 T;

    for (int r=0; r<MAT5_SIZE; ++r)
    {
        for (int c=0; c<MAT5_SIZE; ++c)
        {
            T[r][c] = s * M[r][c];
        }
    }
    
    return T;
}

void eq5 (Matrix5& M, const double s)
{
    for (int r=0; r<MAT5_SIZE; ++r)
    {
        M[r].fill(s);
    }
}

Vector<MAT5_SIZE> mat5Vec5Mul (const Matrix5& M, const Vector<MAT5_SIZE>& V)
{
    Vector<MAT5_SIZE> T;
    T.fill(0.);
    
    for (int r=0; r<MAT5_SIZE; ++r)
    {
        for (int c=0; c<MAT5_SIZE; ++c)
        {
            T[r] += M[r][c] * V[c];
        }
    }
    
    return T;
}

