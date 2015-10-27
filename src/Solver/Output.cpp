#include "Solver.h"

void Solver::log (string fileName)
{
    ofstream out;
    
    //string temps = "log.dat";
    //string slash = "/";
    //string logDir = dir;
    //logDir.append (slash);
    //logDir.append (temps);
    
    out.open (fileName, ofstream::app);
    
    out << endl;
    out << instanceName << endl;
    
    out << "steady = " << steady << endl;
            
    if (!steady)
    {
        out << "dt = " << dt << endl;
        out << "finalTime = " << finalTime << endl;
    }
    
    out << "implicit = " << implicit << endl;
    
    if (implicit)
    {
        if (linearSolverType == 1)
        {
            out << "solver = " << "user GS" << endl;
            out << "nGaussIter = " << nGaussIter << endl;
        }
        else if (linearSolverType == 2)
        {
            out << "solver = PETSC" << endl; // which PETSC? GMRES?
        }
        else
        {
            out << "solver = undefined" << endl; // which PETSC? GMRES?
        }
        
        out << "nImplicitCalls = " << nImplicitCalls << endl;
        out << "cfl = " << cfl << endl;
    }
    
    out << "Temporal order = " << tOrder << endl;
    out << "Spatial order = " << sOrder << endl;
    out << "maxTS/loop = " << maxTimeStep << endl;
    out << "tol = " << tol << endl;
            
    out.close();
}

void Solver::outRes (string fileName)
{
    int s = 10;
    int b = 20;
    
    ofstream out;
    
    if (!out.is_open())
    {
        string dir = fileName;        
        string slash = "/";
        dir.append (slash);
        dir.append (instanceName);
        string temps = "res.dat";        
        dir.append (temps);
        
        out.open (dir, ofstream::app);
    }
        
    if (out.is_open())
    {
        out << left << setw(s) << glo_nTimeStep;
        out << setw(b) << time;
        out << setw(s) << nTimeStep;
        out << setw(b) << aveRes;
        out << endl;
    }
    else
    {
        cout << "could not open file in Solver::outRes(...)" << endl;
        exit(-2);
    }
    
    
    
}