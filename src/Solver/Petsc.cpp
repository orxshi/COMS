#include "Solver.h"

void Solver::Petsc::solveAxb (Grid& gr)
{
    // set values of b
    int ind[vecLocalSize];
    double val[vecLocalSize];
    
    /*int cntr = 0;
    for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
    {
        Cell& cll = gr.cell[ic];
        
        if (cll.iBlank == iBlank_t::FIELD)
        {
            for (int i=0; i<bs; ++i)
            {
                ind[cntr+i] = cntr+i;
                val[cntr+i] = cll.R[i];
            }
            
            cntr += bs;
        }
    }    */
    
    for (PetscInt gp=vecLocBeg; gp<vecLocEnd; ++gp)
    {
        PetscInt brow = static_cast <int> (floor(gp/bs));
        PetscInt c = brow + gr.n_bou_elm;
        PetscInt i = gp % bs;
        
        Cell& cll = gr.cell[c];        
        
        ind[gp-vecLocBeg] = gp;
        val[gp-vecLocBeg] = cll.R[i];
    }
        
    VecSetValues (b, vecLocalSize, ind, val, INSERT_VALUES);
    
    VecAssemblyBegin (b);
    VecAssemblyEnd (b);    
    
    // set values of A
    for (PetscInt brow=matLocBeg/bs; brow<matLocEnd/bs; ++brow)
    {
        PetscInt c = brow + gr.n_bou_elm;
        
        Cell& cll = gr.cell[c];
        PetscInt idxm = brow;
        PetscInt idxn = brow;
        double v[bs*bs];
        
        for (int q=0; q<bs; ++q)
        {
            for (int w=0; w<bs; ++w)
            {
                v[q*bs+w] = cll.D[q][w];
            }
        }
        
        MatSetValuesBlocked (A, 1, &idxm, 1, &idxn, v, INSERT_VALUES);        
        
        for (int nn=0; nn<cll.nei.size(); ++nn)
        {
            if (cll.nei[nn] >= gr.n_bou_elm)
            {
                idxn = cll.nei[nn] - gr.n_bou_elm;
                
                Face& f = gr.face[cll.face[nn]];
                if (c == f.nei[0])
                {
                    for (int q=0; q<bs; ++q)
                    {
                        for (int w=0; w<bs; ++w)
                        {
                            v[q*bs+w] = f.M[1][q][w];
                        }
                    }
                }
                else
                {
                    for (int q=0; q<bs; ++q)
                    {
                        for (int w=0; w<bs; ++w)
                        {
                            v[q*bs+w] = -f.M[0][q][w];
                        }
                    }
                }
                
                MatSetValuesBlocked (A, 1, &idxm, 1, &idxn, v, INSERT_VALUES);
            }
        }
    }
    
    /*// set values of A
    for (PetscInt gp=first; gp<last; ++gp)
    {
        PetscInt brow = static_cast <int> (floor(gp/bs));
        PetscInt gpDiagStart = brow*bs;
        PetscInt c = brow + gr.n_bou_elm;
        PetscInt i = gp % bs;
        
        Cell& cll = gr.cell[c];
        PetscInt idxn[bs];
        double v[bs];
        
        for (int j=0; j<bs; ++j) { idxn[j] = gpDiagStart + j; }        
        for (int j=0; j<bs; ++j) { v[j] = cll.D[i][j]; }
        
        MatSetValues (A, 1, &gp, bs, idxn, v, INSERT_VALUES);
        
        for (int nn=0; nn<cll.nei.size(); ++nn)
        {
            if (cll.nei[nn] >= gr.n_bou_elm)
            {
                PetscInt neiBcol = cll.nei[nn] - gr.n_bou_elm;
                PetscInt gpNeiStart = neiBcol * bs;
                
                for (int j=0; j<bs; ++j) { idxn[j] = gpNeiStart + j; }
                
                Face& f = gr.face[cll.face[nn]];
                if (c == f.nei[0])
                {
                    for (int j=0; j<bs; ++j) { v[j] = f.M[1][i][j]; }
                }
                else
                {
                    for (int j=0; j<bs; ++j) { v[j] = -f.M[0][i][j]; }
                }
                
                MatSetValues (A, 1, &gp, bs, idxn, v, INSERT_VALUES);
            }
        }
    }*/
    
    MatAssemblyBegin (A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY);    
    //
    
    KSPSolve (ksp, b, x);
    
    //
    int recvcounts[nProcs];
    int displs[nProcs];
    displs[0] = 0;
    recvcounts[0] = localSizes[0];
    
    for (int i=1; i<nProcs; ++i)
    {
        recvcounts[i] = localSizes[i];
        displs[i] = displs[i-1] + recvcounts[i-1];
    }
    
    double *dx = NULL;
    VecGetArray (x, &dx);
    MPI_Allgatherv (dx, localSizes[rank], MPI_DOUBLE, DX, localSizes, displs, MPI_DOUBLE, world);
    VecRestoreArray (x, &dx);
    //
    
    for (PetscInt gp=0; gp<vecGlobalSize; ++gp)
    {
        PetscInt brow = static_cast <int> (floor(gp/bs));
        PetscInt c = brow + gr.n_bou_elm;
        PetscInt i = gp % bs;        
        Cell& cll = gr.cell[c];
        
        cll.dQ[i] = DX[gp];
    }
}

void Solver::Petsc::finalize()
{
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    //PetscFree (DX);
    delete DX;
    //free (DX);
    DX = NULL;
    delete [] localSizes;
}
