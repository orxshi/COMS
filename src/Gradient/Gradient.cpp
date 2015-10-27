#include "../Grid/Grid.h"

void Grid::leastSquaresCoeffs()
{
    double dx, dy, dz;
    CVector d;

    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        Cell& cll = cell[c];
        
        for (const int f: cell[c].nei)
        //for (const Cell& f: cll.nei)
        {            
            d = cell[f].cnt - cll.cnt;

            dx = d[0];
            dy = d[1];
            dz = d[2];

            cll.r_11 += dx * dx;
            cll.r_12 += dx * dy;
            cll.r_13 += dx * dz;
            cll.r_22 += dy * dy;
            cll.r_23 += dy * dz;
            cll.r_33 += dz * dz;
        }

        cll.r_11 = sqrt (cll.r_11);
        cll.r_12 /= cll.r_11;
        cll.r_13 /= cll.r_11;
        cll.r_22 = sqrt (cll.r_22 - pow(cll.r_12,2));
        cll.r_23 = (cll.r_23 - cll.r_12 * cll.r_13) / cll.r_22;
        cll.r_33 = sqrt (cll.r_33 - pow(cll.r_13,2) - pow(cll.r_23,2));
    }
}

void Grid::leastSquaresGrad()
{
    double dx, dy, dz, a1, a2, a3, psi, Wx, Wy, Wz;
    double tempf;
    CVector d;

    for (int c=n_bou_elm; c<cell.size(); ++c)
    //for (Cell& cll: cell)
    {
        Cell& cll = cell[c];
        
        for (int i=0; i<N_VAR; ++i)
        {
            cll.grad[i].fill(0.);
        }
        
        for (const int f: cll.nei)
        {
            d = cell[f].cnt - cll.cnt;

            dx = d[0];
            dy = d[1];
            dz = d[2];

            a1 = dx / pow (cll.r_11,2);
            a2 = (dy - dx * cll.r_12 / cll.r_11) / pow (cll.r_22,2);
            psi = (cll.r_12 * cll.r_23 - cll.r_13 * cll.r_22) / (cll.r_11 * cll.r_22);
            a3 = (dz - dy * cll.r_23 / cll.r_22 + psi * dx) / pow (cll.r_33,2);
            Wx = a1 - a2 * cll.r_12 / cll.r_11 + psi * a3;
            Wy = a2 - a3 * cll.r_23 / cll.r_22;
            Wz = a3;

            for (int n=0; n<N_VAR; ++n)
            {
                tempf = cell[f].prim[n] - cll.prim[n];

                cll.grad[n][0] += Wx * tempf;
                cll.grad[n][1] += Wy * tempf;
                cll.grad[n][2] += Wz * tempf;
            }
        }
    }
}