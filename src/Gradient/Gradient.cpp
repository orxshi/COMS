#include "Gradient.h"

Gradient::Gradient (Grid& gr)
{
    grad.resize (gr.n_in_elm);
    
    for (int ic=0; ic<gr.n_in_elm; ++ic)
    {
        for (int i=0; i<N_VAR; ++i)
        {
            grad[ic][i].fill(0.);
        }
    }
    
    Wx.resize (gr.n_in_elm);
    Wy.resize (gr.n_in_elm);
    Wz.resize (gr.n_in_elm);

    leastSquaresCoeffs (gr);
    initParallelVars (gr);
}

void Gradient::initParallelVars (Grid& gr)
{
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    localSizes  = new int [nProcs];
    localSizesNVARNDIM  = new int [nProcs];
    displs = new int [nProcs];
    displsNVARNDIM = new int [nProcs];

    int rem;

    // get localSizeFace
    rem = gr.n_in_elm % nProcs;
    
    //if (rem != 0)
    //{        
        localSize = (gr.n_in_elm - rem) / nProcs;
    //}
    
    if (rank == nProcs-1)
    {
        localSize += rem;
    }
    
    localSizeNVARNDIM = localSize * N_VAR * N_DIM;
    
    // gather localSizesFace and localSizesFaceM, F, V
    MPI_Allgather (&localSize, 1, MPI_INT, localSizes, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather (&localSizeNVARNDIM, 1, MPI_INT, localSizesNVARNDIM, 1, MPI_INT, MPI_COMM_WORLD);
    
    displs[0] = gr.n_bou_elm;
    for (int i=1; i<nProcs; ++i)
    {
        displs[i] = displs[i-1] + localSizes[i-1];
    }
    
    displsNVARNDIM[0] = 0;
    for (int i=1; i<nProcs; ++i)
    {
        displsNVARNDIM[i] = displsNVARNDIM[i-1] + localSizesNVARNDIM[i-1];
    }
}

void Gradient::leastSquaresCoeffs (Grid& gr)
{    
    double dx, dy, dz, a1, a2, a3, psi;
    double r_11;
    double r_12;
    double r_13;
    double r_22;
    double r_23;
    double r_33;
    CVector d;

    for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
    {
        Cell& cll = gr.cell[c];
        
        for (const int f: gr.cell[c].nei)        
        {            
            d = gr.cell[f].cnt - cll.cnt;

            dx = d[0];
            dy = d[1];
            dz = d[2];

            r_11 += dx * dx;
            r_12 += dx * dy;
            r_13 += dx * dz;
            r_22 += dy * dy;
            r_23 += dy * dz;
            r_33 += dz * dz;
        }

        r_11 = sqrt (r_11);
        r_12 /= r_11;
        r_13 /= r_11;
        r_22 = sqrt (r_22 - pow(r_12,2));
        r_23 = (r_23 - r_12 * r_13) / r_22;
        r_33 = sqrt (r_33 - pow(r_13,2.) - pow(r_23,2.));
        
        for (const int f: gr.cell[c].nei)
        {
            d = gr.cell[f].cnt - cll.cnt;

            dx = d[0];
            dy = d[1];
            dz = d[2];

            a1 = dx / pow (r_11,2.);
            a2 = (dy - dx * r_12 / r_11) / pow (r_22,2.);
            psi = (r_12 * r_23 - r_13 * r_22) / (r_11 * r_22);
            a3 = (dz - dy * r_23 / r_22 + psi * dx) / pow (r_33,2.);
            
            Wx[c-gr.n_bou_elm].push_back (a1 - a2 * r_12 / r_11 + psi * a3);
            Wy[c-gr.n_bou_elm].push_back (a2 - a3 * r_23 / r_22);
            Wz[c-gr.n_bou_elm].push_back (a3);
        }
    }    
}

void Gradient::leastSquaresGrad (Grid& gr)
{    
    double tempf;
    CVector d;

    for (int ic=displs[rank]; ic<displs[rank]+localSize; ++ic)
    //for (int c=n_bou_elm; c<cell.size(); ++c)
    //for (Cell& cll: cell)
    {
        Cell& cll = gr.cell[ic];
        
        for (int i=0; i<N_VAR; ++i)
        {
            //cll.grad[i].fill(0.);
            grad[ic-gr.n_bou_elm][i].fill(0.);
        }
        
        for (const int f: cll.nei)
        {
            for (int n=0; n<N_VAR; ++n)
            {
                tempf = gr.cell[f].prim[n] - cll.prim[n];

                grad[ic-gr.n_bou_elm][n][0] += Wx[ic-gr.n_bou_elm][n] * tempf;
                grad[ic-gr.n_bou_elm][n][1] += Wy[ic-gr.n_bou_elm][n] * tempf;
                grad[ic-gr.n_bou_elm][n][2] += Wz[ic-gr.n_bou_elm][n] * tempf;
            }
        }
    }
    
    MPI_Allgatherv (MPI_IN_PLACE, localSizesNVARNDIM[rank], MPI_DOUBLE, &grad[0][0][0], localSizesNVARNDIM, displsNVARNDIM, MPI_DOUBLE, MPI_COMM_WORLD);
}