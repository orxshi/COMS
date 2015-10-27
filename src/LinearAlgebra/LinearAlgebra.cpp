#include "LinearAlgebra.h"

void GaussElim (Vector2D<N_DIM,N_DIM> A, CVector b, CVector& x)
{
    /****************************************************

        Solve linear system Ax = b
        using Gaussian elimination without pivoting
        A is an n by n matrix
        b is an n by k matrix (k copies of n-vectors)
        x is an n by k matrix (k copies of solution vectors)

        http://cis.poly.edu/~mleung/CS3734/s03/ch02/Gauss.htm

    ****************************************************/

    double m, tempB;

    for (unsigned int i=0; i<x.size()-1; ++i)
    {
        for (unsigned int j=i+1; j<x.size(); ++j)
        {
            m = - A[j][i] / A[i][i];

            for(unsigned int k=i; k<x.size(); ++k)
            {
                A[j][k] += (m * A[i][k]);
            }

            b[j] += m * b[i];
        }
    }

    x[x.size() - 1] = b[x.size() - 1] / A[x.size() - 1][x.size() - 1];

    // Use back substitution to find unknowns
    for (int i=x.size()-2; i>=0; --i)
    {
        tempB = b[i];

        for(unsigned j=i+1; j<x.size(); ++j)
        {
            tempB -= (A[i][j] * x[j]);
        }

        x[i] = tempB / A[i][i];
    }
}

void newtonSolve(const Vector2D<N_DIM,8>& f, double& u, double& v, double& w)
{
    /**********************************************************

        Solve the non-linear vector equation
        f(u,v,w)=0
        using Newton-Raphson Method

        f(u,v,w)=f0+f1*u+f2*v+f3*w+f4*u*v+f5*u*w+f6*v*w+f7*u*v*w

        Algorithm:

        q=[u,v,w]
        iterate:
        dq=-[df/dq]^-1*f(q)

    **********************************************************/

    int itmax = 500;
    double uv, wu, vw, uvw, norm, convergenceLimit = 1e-14;
    Vector<3> rhs, x;
    Vector2D<3,3> lhs;
    double alph = 1.;

    u = 0.5;
    v = 0.5;
    w = 0.5;

    for (int i=0; i<itmax; ++i)
    {
        uv  = u * v;
        vw  = v * w;
        wu  = w * u;
        uvw = u * v * w;

        for (int j=0; j<3; ++j)
        {
            rhs[j] = f[0][j] +
                     f[1][j] * u +
                     f[2][j] * v +
                     f[3][j] * w +
                     f[4][j] * uv +
                     f[5][j] * vw +
                     f[6][j] * wu +
                     f[7][j] * uvw;
        }
        
        norm = pow(rhs[0],2) + pow(rhs[1],2) + pow(rhs[2],2);

        if (sqrt(norm) <= convergenceLimit)
        {
            break;
        }

        for (int j=0; j<3; ++j)
        {
            lhs[j][0] = f[1][j] + f[4][j] * v + f[6][j] * w + f[7][j] * vw;
            lhs[j][1] = f[2][j] + f[5][j] * w + f[4][j] * u + f[7][j] * wu;
            lhs[j][2] = f[3][j] + f[6][j] * u + f[5][j] * v + f[7][j] * uv;
        }

        GaussElim (lhs, rhs, x);

        u = u - x[0] * alph;
        v = v - x[1] * alph;
        w = w - x[2] * alph;
    }
}

void osInterpolants (const vector<CVector>& xc, const CVector& xp, const unsigned int nVtx, double frac[])
{
    //--declarations----------------------------
    double u, v, w;
    double oneminusU;
    double oneminusV;
    double oneminusW;
    double oneminusUV;
    Vector<3> rhs, x;
    Vector2D<3,8> f;
    Vector2D<3,3> lhs;
    //------------------------------------------

    if (nVtx == 4)
    {
        for (unsigned int i=0; i<3; ++i)
        {
            rhs[i] = xp[i] - xc[3][i];
        }

        for (unsigned int i=0; i<3; ++i)
        {
            for (unsigned int j=0; j<3; ++j)
            {
                lhs[j][i] = xc[j][i] - xc[4][i];
            }

        }

        GaussElim (lhs, rhs, x);

        frac[0] = x[0];
        frac[1] = x[1];
        frac[2] = x[2];
        frac[3] = 1. - (frac[0] + frac[1] + frac[2]);
    }
    else if (nVtx == 6) // prism
    {
        for (unsigned int i=0; i<3; ++i)
        {
            f[0][i] = xc[0][i] - xp[i];
            f[1][i] = xc[1][i] - xc[0][i];
            f[2][i] = xc[2][i] - xc[0][i];
            f[3][i] = xc[3][i] - xc[0][i];
            f[4][i] = 0.;
            f[5][i] = xc[0][i] - xc[2][i] - xc[3][i] + xc[5][i];
            f[6][i] = xc[0][i] - xc[1][i] - xc[3][i] + xc[4][i];
            f[7][i] = 0.;
        }

        newtonSolve (f, u, v, w);

        oneminusUV = 1. - u - v;
        oneminusU  = 1. - u;
        oneminusV  = 1. - v;
        oneminusW  = 1. - w;

        frac[0] = oneminusUV*oneminusW;
        frac[1] = u*oneminusW;
        frac[2] = v*oneminusW;
        frac[3] = oneminusUV*w;
        frac[4] = u*w;
        frac[5] = v*w;
    }
    else if (nVtx == 8) // hexahedron
    {
        for (unsigned int i=0; i<3; ++i)
        {
            f[0][i] = xc[0][i] - xp[i];
            f[1][i] = xc[1][i] - xc[0][i];
            f[2][i] = xc[3][i] - xc[0][i];
            f[3][i] = xc[4][i] - xc[0][i];
            f[4][i] = xc[0][i] - xc[1][i] + xc[2][i] - xc[3][i];
            f[5][i] = xc[0][i] - xc[3][i] + xc[7][i] - xc[4][i];
            f[6][i] = xc[0][i] - xc[1][i] + xc[5][i] - xc[4][i];
            f[7][i] = -xc[0][i] + xc[1][i] - xc[2][i] + xc[3][i] + xc[4][i] - xc[5][i] + xc[6][i] - xc[7][i];
        }

        newtonSolve (f, u, v, w);

        oneminusU = 1. - u;
        oneminusV = 1. - v;
        oneminusW = 1. - w;

        frac[0] = oneminusU*oneminusV*oneminusW;
        frac[1] = u*oneminusV*oneminusW;
        frac[2] = u*v*oneminusW;
        frac[3] = oneminusU*v*oneminusW;
        frac[4] = oneminusU*oneminusV*w;
        frac[5] = u*oneminusV*w;
        frac[6] = u*v*w;
        frac[7] = oneminusU*v*w;
    }
}
