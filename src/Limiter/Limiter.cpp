#include "Limiter.h"

Limiter::Limiter (Grid& gr)
{
    // defaults
    type = 3;
    
    if (type == 2)
    {
        ksiBJ.resize(gr.n_in_elm);
        for (int i=0; i<ksiBJ.size(); ++i)
        {
            ksiBJ[i].resize(N_VAR);
        }
    }
    
    if (type == 3)
    {
        ksiV.resize(gr.n_in_elm);
        for (int i=0; i<ksiV.size(); ++i)
        {
            ksiV[i].resize(N_VAR);
        }
    }
    
    initParallelVars (gr);
}

void Limiter::initParallelVars (Grid& gr)
{
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    localSizes  = new int [nProcs];
    localSizesNVAR  = new int [nProcs];
    displs = new int [nProcs];
    displsNVAR = new int [nProcs];

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
    
    localSizeNVAR = localSize * N_VAR;
    
    // gather localSizesFace and localSizesFaceM, F, V
    MPI_Allgather (&localSize, 1, MPI_INT, localSizes, 1, MPI_INT, MPI_COMM_WORLD);    
    MPI_Allgather (&localSizeNVAR, 1, MPI_INT, localSizesNVAR, 1, MPI_INT, MPI_COMM_WORLD);    
    
    displs[0] = gr.n_bou_elm;
    for (int i=1; i<nProcs; ++i)
    {
        displs[i] = displs[i-1] + localSizes[i-1];
    }
    
    displsNVAR[0] = 0;
    for (int i=1; i<nProcs; ++i)
    {
        displsNVAR[i] = displsNVAR[i-1] + localSizesNVAR[i-1];
    }
}

void Limiter::getLimitedGrad (const Vector2D<3,N_VAR>& gradL, const Vector2D<3,N_VAR>& gradR, Vector2D<3,N_VAR>& grad)
{
    double R;
    double phi = 0.;

    for (int k=0; k<N_VAR; ++k)
    {
        for (int d=0; d<N_DIM; ++d)
        {
            if (gradL[k][d] != 0.)
            {
                R = gradR[k][d] / gradL[k][d];
            }
            else
            {
                R = -1.;
            }

            if ( R >= 0. )
            {
                if (type == 0)
                {
                    double tmp = 2. / (1. + R);            
                    phi = min ( tmp, R*tmp );
                }
                else if (type == 1)
                {
                    phi = (2. * R) / (pow(R,2.) + 1.);
                }
                
                grad[k][d] = 0.5 * phi * (gradL[k][d] + gradR[k][d]);
            }
            else
            {
                grad[k][d] = 0.;
            }
        }
    }
}

void Limiter::getLimitedGradDarwish (const Vector2D<3,N_VAR>& gradL, const Vector2D<3,N_VAR>& gradR, const CVector& rL, const CVector& rR,
                                const Vector<N_VAR>& varL, const Vector<N_VAR>& varR, Vector<N_VAR>& reconstL, Vector<N_VAR>& reconstR)
{
    double rf;
    double ksi;
    Vector2D<3,N_VAR> gradC;
    Vector<N_VAR> varC;
    Vector<N_VAR> varD;
    Vector<N_VAR> reconstC;
    Vector<N_VAR> reconstD;
    CVector rCD;
    
    CVector rLR = rR - rL;
    CVector vLR;
    vLR[0] = varR[0] - varL[0];
    vLR[1] = varR[1] - varL[1];
    vLR[2] = varR[2] - varL[2];
    double dire = dotP(vLR,rLR);
    
    if (dire > 0.)
    {
        // L is upwind
        // R is downwind        
        gradC = gradL;
        varC = varL;
        varD = varR;
        rCD = rLR;
    }
    else if (dire < 0.)
    {
        // R is upwind
        // L is downwind        
        gradC = gradR;
        varC = varR;
        varD = varL;
        rCD = -1. * rLR;
    }
    else
    {
        // grad is zero, no need for limiter
        reconstL = varL;
        reconstR = varR;
        return;
    }

    for (int k=0; k<N_VAR; ++k)
    {
        rf = 2. * dotP(gradC[k],rCD) / (varD[k] - varC[k]) - 1.;
        
        if (type == 0)
        {
            ksi = max (0., min(1.,rf));
        }
        
        reconstC[k] = varC[k] + 0.5 * ksi * (varD[k] - varC[k]);
        reconstD[k] = varD[k] - 0.5 * ksi * (varD[k] - varC[k]);
    }
    
    if (dire > 0.)
    {
        reconstL = reconstC;
        reconstR = reconstD;
    }
    else if (dire < 0.)
    {
        reconstL = reconstD;
        reconstR = reconstC;
    }
}

void Limiter::bj (Grid& gr)
{    
    /*Vector<N_VAR> ksiF;
    int icc;

    for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
    {
        Cell& cll = gr.cell[ic];
        icc = ic - gr.n_bou_elm;
        
        for (int k=0; k<N_VAR; ++k)
        {
            ksiBJ[icc][k] = BIG_POS_NUM;
        }
    
        for (int iFace: cll.face)
        {
            Face& f = gr.face[iFace];
            
            int iNei;
            for (int i=0; i<f.nei.size(); ++i)
            {                
                if (f.nei[i] != ic)
                {
                    iNei = f.nei[i];
                    break;
                }
            }            
            Cell& nei = gr.cell[iNei];
            
            CVector dis = f.cnt - cll.cnt;
            
            for (int k=0; k<N_VAR; ++k)
            {
                double tmp = dotP(cll.grad[k],dis);
            
                if (tmp > 0.)
                {
                    ksiF[k] = min( (max( nei.prim[k], cll.prim[k] ) - cll.prim[k]) / tmp, 1. );
                }
                else if (tmp < 0.)
                {
                    ksiF[k] = min( (min( nei.prim[k], cll.prim[k] ) - cll.prim[k]) / tmp, 1. );
                }
                else
                {
                    ksiF[k] = 1.;
                }
                
                ksiBJ[icc][k] = min (ksiBJ[icc][k], ksiF[k]);
            }
        }
    }*/
}

void Limiter::venka (Grid& gr, Gradient& gradient)
{    
    auto phi = [&] (double x, double eps)
    {
        return (pow(x,2.) + 2.*x + eps) / (pow(x,2.) + x + 2. + eps);
    };

    Vector<N_VAR> ksiF;
    int icc;
    double K = 0.3;
    
    //for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
    for (int ic=displs[rank]; ic<displs[rank]+localSize; ++ic)
    {
        Cell& cll = gr.cell[ic];
        icc = ic - gr.n_bou_elm;
        double eps = pow(K*cll.vol, 3.);
        
        for (int k=0; k<N_VAR; ++k)
        {
            ksiV[icc][k] = BIG_POS_NUM;
        }
    
        for (int iFace: cll.face)
        {
            Face& f = gr.face[iFace];
            
            int iNei;
            for (int i=0; i<f.nei.size(); ++i)
            {                
                if (f.nei[i] != ic)
                {
                    iNei = f.nei[i];
                    break;
                }
            }            
            Cell& nei = gr.cell[iNei];
            
            CVector dis = f.cnt - cll.cnt;
            
            for (int k=0; k<N_VAR; ++k)
            {                
                //double tmp = dotP(cll.grad[k],dis);
                double tmp = dotP(gradient.grad[icc][k],dis);                
            
                if (tmp > 0.)
                {                    
                    ksiF[k] = phi ((max( nei.prim[k], cll.prim[k] ) - cll.prim[k]) / tmp, eps);
                }
                else if (tmp < 0.)
                {
                    ksiF[k] = phi ((min( nei.prim[k], cll.prim[k] ) - cll.prim[k]) / tmp, eps);
                }
                else
                {
                    ksiF[k] = 1.;
                }
                
                ksiV[icc][k] = min (ksiV[icc][k], ksiF[k]);
            }
        }
    }
    
    MPI_Allgatherv (MPI_IN_PLACE, localSizesNVAR[rank], MPI_DOUBLE, &ksiV[0][0], localSizesNVAR, displsNVAR, MPI_DOUBLE, MPI_COMM_WORLD);
}

void minMod(const Vector2D<3,N_VAR>& gradL, const Vector2D<3,N_VAR>& gradR, Vector2D<3,N_VAR>& grad)
{
    // http://www.cfdbooks.com/cfdcodes/oned_euler_v1.f90

    double tmp;

    for (int k=0; k<N_VAR; ++k)
    {
        for (int d=0; d<N_DIM; ++d)
        {
            tmp = gradL[k][d] * gradR[k][d];

            if ( tmp > 0. )
            {
                if ( fabs(gradL[k][d]) < fabs(gradR[k][d]) )
                {
                    grad[k][d] = gradL[k][d];
                }
                else
                {
                    grad[k][d] = gradR[k][d];
                }
            }
            else
            {
                grad[k][d] = 0.;
            }
        }
    }
}





double venkata (const Vector2D<3,N_VAR>& grad, Vector<N_VAR>& uMax, Vector<N_VAR>& uMin, Vector<N_VAR>& u, const CVector& dis)
{
    double ksi;
    Vector<N_VAR> d;
    
    auto phi = [&] (double y) { return (pow(y,2.) + 2.*y) / (pow(y,2.) + y + 2.); };
    
    for (int i=0; i<N_VAR; ++i)
    {
        d[i] = dotP (grad[i], dis);
        
        double argMax = (uMax[i] - u[i]) / d[i];
        double argMin = (uMin[i] - u[i]) / d[i];
        
        if (d[i] > 0.)
        {
            ksi = phi (argMax);
        }
        else if (d[i] < 0.)
        {
            ksi = phi (argMin);
        }
        else
        {
            ksi = 1.;
        }
    }
    
    return ksi;
}
