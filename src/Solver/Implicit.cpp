#include "Solver.h"

void Solver::impl (Grid& gr)
{    
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    //Watch wt;
    
    //preSolverCheck (gr);
    
    string dir = gr.outputDir;
    string temps = "res.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);
    int printThres = 100;
    
    Limiter limiter (gr);
    
    if (sOrder == 2)
    {
        //gr.leastSquaresGrad();
        gradient.leastSquaresGrad (gr);
    }
    roe.roeflx (gr, limiter, M0, M1, gradient); // parallel
    
    for (nTimeStep=0; nTimeStep<maxTimeStep; ++nTimeStep)
    {        
        interflux(gr); // serial
        
        switch (linearSolverType)
        {
            case 1:
                gauss_seidel (gr); // only fields
                break;
            case 2:                
                petsc.solveAxb (gr, M0, M1);                
                break;
            default:
                cout << "undefined linear solver in Solver::impl(...)" << endl;
                exit(-2);
                break;
        }
        
        diff_to_cons_prim (gr); // only fields
        getRes (gr, limiter); // only fields
        
        if (rank == MASTER_RANK) { outRes(gr.outputDir); }
        gr.apply_BCs();
        
        if (verbose && rank == MASTER_RANK && nTimeStep%printThres==0)
        {
            cout << left << setw(10) << fixed << time;
            cout << setw(10) << nTimeStep;
            cout << scientific << aveRes << endl;            
        }        
        
        if (fabs(aveRes) < tol) { break; }
        
        ++glo_nTimeStep;
    }
    
    for (Cell& cll: gr.cell)
    {
        cll.oldold_cons = cll.old_cons;
        cll.old_cons = cll.cons;
    }    
    
    ++nImplicitCalls;
}
