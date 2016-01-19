/* 
 * File:   Solver.h
 * Author: Orhan Shibliyev
 *
 * Created on October 2, 2014, 10:00 PM
 */

#ifndef SOLVER_H
#define	SOLVER_H

//#include <petscsys.h>
#include <petscksp.h>
//#include "../Grid/Grid.h"
#include "../Output/Output.h"
#include "../InterFlux/Roe/Roe.h"
#include "../Time/Time.h"
#include "../Gradient/Gradient.h"

//enum class tsOrder_t {UNDEFINED=0, FIRST=1, SECOND=2}; // temporal and spatial orders
//enum class LinearSol_t {UNDEFINED=0, MYGS=1, PETSC=2}; // linear solver for implicit method

struct Solver
{
    struct Petsc
    {
        Vec x, b;
        Mat A;
        KSP ksp;
        PC pc;        
        PetscInt n;
        PetscInt bs;
        PetscInt vecGlobalSize;
        PetscInt vecLocalSize;
        PetscInt vecLocBeg, vecLocEnd;
        PetscInt matLocBeg, matLocEnd;
        double* DX;                
        int* localSizes;
        int matLocalSize;
        
        Petsc (Grid& gr);
        void solveAxb (Grid& gr, vector <Matrixd<N_VAR,N_VAR>>& M0, vector <Matrixd<N_VAR,N_VAR>>& M1);
        void finalize();
    } petsc;    
    
    int nGaussIter;
    int maxTimeStep;
    int nTimeStep;
    int glo_nTimeStep;
    int nImplicitCalls;
    int tOrder;
    int sOrder;
    int linearSolverType;
    int nSampleImpl;
    int nSampleInne;    
    double cfl;
    double dt;
    double finalTime;
    double tol;    
    double time;
    double aveRes;
    double eTimeImpl;
    double eTimeInne;
    bool steady;
    bool implicit;
    bool verbose;
    bool isSampledImpl;
    bool isSampledInne;
    //tsOrder_t tOrder;
    //tsOrder_t sOrder;
    //LinearSol_t linearSolverType;        
    string instanceName;
    //Matrixd<N_VAR, N_VAR> R;
    //Matrixd<N_VAR, N_VAR> D;
    //MPI_Comm world;
    Gradient gradient;
    Limiter limiter;
    Watch wImpl;
    Watch wInne;
    
    vector <Matrixd<N_VAR,N_VAR>> M0;
    vector <Matrixd<N_VAR,N_VAR>> M1;
    Roe roe;
    
    Solver (Grid& gr, string instanceName);
    
    void expl (Grid& gr); // don't use
    void impl (Grid& gr);
    void interflux (Grid& gr);
    void gauss_seidel (Grid& g);
    void updateVars (Grid& gr);
    double setExpRes (Grid& gr); // dont use
    //void preSolver(Grid& gr);
    void log (string fileName);
    void outRes (string fileName);
    void preSolverCheck(const Grid& gr);    
    void read (string fileName);
    void set_residual (Grid& g); // don't use
    void diff_to_cons_prim(Grid& g);
    void getRes (Grid& gr);    
    void printTimeInne();
    void printTimeImpl();
    
    bool cm (string s, ifstream& in);
    template <class T> void cmh (string s, string membs, T& memb, ifstream& in, bool& found) // helper fnc for cm
    {
        if (!found)
        {
            int isEqual = 0;

            if ( s.compare ( membs ) == isEqual )
            {
                in >> memb;
                found = true;
                
                /*cout << s << endl;
                cout << membs << endl;  
                cout << memb << endl;*/
            }
        }
    }
};

#endif	/* SOLVER_H */

