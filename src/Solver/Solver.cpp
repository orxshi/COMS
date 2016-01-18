#include "Solver.h"

Solver::Solver (Grid& gr, string instanceName) : petsc(gr), roe(gr), gradient (gr), limiter (gr)
{
    //default    
    tOrder = 2;
    sOrder = 1;
    linearSolverType = 1; // MYGS
    nGaussIter = 5;
    maxTimeStep = 10000;
    cfl = 5;
    dt = 1.;
    finalTime = 10.;
    tol = 1e-12;    
    steady = false;
    implicit = true;
    verbose = true;
    
    // initialize
    time = 0.;
    glo_nTimeStep = 0;
    nImplicitCalls = 0;    
    this->instanceName = instanceName;
    M0.resize (gr.face.size());
    M1.resize (gr.face.size());
}



Solver::Petsc::Petsc (Grid& gr)
{    
    int nProcs;

    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    /*n = 0;
    for (Cell& cll: gr.cell)
    {
        if (cll.iBlank == iBlank_t::FIELD)
        {
            ++n;
        }
    }*/
    
    n = gr.n_in_elm;
    bs = N_VAR;
    vecGlobalSize = n*bs;
    //DX = (double*) malloc (vecGlobalSize*sizeof(double));
    //PetscMalloc1 (xGlobalSize, &DX);
    DX = new double [vecGlobalSize];
    
    
    
    // set x
    VecCreate (world, &x);
    VecSetType (x, VECSTANDARD);
    VecSetBlockSize(x, bs);
    VecSetSizes (x, PETSC_DECIDE, vecGlobalSize);
    VecGetLocalSize (x, &vecLocalSize);
    VecGetOwnershipRange (x, &vecLocBeg, &vecLocEnd);
    
    

    // set b
    VecDuplicate (x, &b);
    
    // set A
    MatCreate (world, &A);    
    MatSetType (A, MATMPIBAIJ);
    //MatSetType (A, MATMPIAIJ);
    MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE, vecGlobalSize, vecGlobalSize);
    //MatSeqBAIJSetPreallocation (A, bs, 4, NULL);
    //MatSeqAIJSetPreallocation (A, 4*bs, NULL);    
    //MatMPIAIJSetPreallocation (A, 4*bs, NULL, 4*bs, NULL);
    MatMPIBAIJSetPreallocation (A, bs, 5, NULL, 5, NULL);
    //MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    //MatCreateBAIJ (world, bs, PETSC_DECIDE, PETSC_DECIDE, n*bs, n*bs, 1, NULL, 3, NULL, &A);
    // specific to pentagonal mesh . change later
    MatGetOwnershipRange (A, &matLocBeg, &matLocEnd);
    MatGetLocalSize (A, &matLocalSize, NULL);
    
    
    
    KSPCreate (world, &ksp);
    KSPSetOperators (ksp, A, A);
    KSPGetPC (ksp, &pc);
    PCSetType (pc, PCSOR);
    KSPSetType (ksp, KSPGMRES);
    KSPSetFromOptions (ksp);
        
    localSizes = new int [nProcs];
    int recvcounts[nProcs];
    int displs[nProcs];
    displs[0] = 0;    
    for (int i=0; i<nProcs; ++i)
    {
        recvcounts[i] = 1;
    }
    for (int i=1; i<nProcs; ++i)
    {
        displs[i] = displs[i-1] + recvcounts[i-1];
    }
    
    MPI_Allgatherv (&vecLocalSize, 1, MPI_INT, localSizes, recvcounts, displs, MPI_INT, world);
}

void Solver::preSolverCheck (const Grid& gr)
{
    //string reason = "nothing";
    
    for (const Cell& cll: gr.cell)
    {
        if (cll.iBlank == iBlank_t::UNDEFINED)
        {
            cout << "iBlank is undefined in Solver::preSolverCheck(...)" << endl;
            exit(-2);
        }
    }
}

bool Solver::cm (string s, ifstream& in)
{    
    bool found = false;
    
    cmh (s, STRINGTIFY(nGaussIter), nGaussIter, in, found);
    cmh (s, STRINGTIFY(maxTimeStep), maxTimeStep, in, found);
    cmh (s, STRINGTIFY(cfl), cfl, in, found);
    cmh (s, STRINGTIFY(dt), dt, in, found);
    cmh (s, STRINGTIFY(finalTime), finalTime, in, found);
    cmh (s, STRINGTIFY(tol), tol, in, found);    
    cmh (s, STRINGTIFY(steady), steady, in, found);
    cmh (s, STRINGTIFY(implicit), implicit, in, found);
    cmh (s, STRINGTIFY(verbose), verbose, in, found);
    cmh (s, STRINGTIFY(tOrder), tOrder, in, found);
    cmh (s, STRINGTIFY(sOrder), sOrder, in, found);
    cmh (s, STRINGTIFY(linearSolverType), linearSolverType, in, found);
    
    return found;
}

void Solver::read (string fileName)
{
    string tmps;
    ifstream in;
    in.open (fileName);
    
    if (in.is_open())
    {
        in >> tmps;
        while ( !in.eof() )
        {
            if ( cm (tmps, in) == false )
            {
                cout << "undefined input in Solver::read(...)" << endl;
                exit(-2);
            }

            in >> tmps;
        }
    }
    else
    {
        cout << "could not open file in Solver::read(...)" << endl;
        exit(-2);
    }    

    in.close();
}