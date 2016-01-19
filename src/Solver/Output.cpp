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
    
    if (nImplicitCalls < nSampleImpl)
    {
        out << "could not get enough number of samples for implicit solver" << endl;
        out << "expected number of samples = " << nSampleImpl << endl;
        out << "real number of samples = " << nImplicitCalls << endl;
    }
    else
    {
        cout << "time spent for implicit solver for " << nSampleImpl << " samples = " << fixed << setprecision(2) << (eTimeImpl /= nSampleImpl) << wImpl.unit << endl;
    }
    
    cout << "time spent for inner iter for " << nSampleInne << " samples = " << fixed << setprecision(2) << (eTimeInne /= nSampleInne) << wInne.unit << endl;
            
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

void Solver::printTimeInne()
{    
    cout << "time spent for inner iter for " << nSampleInne << " samples = " << fixed << setprecision(2) << (eTimeInne /= nSampleInne) << wInne.unit << endl;    
}

void Solver::printTimeImpl()
{   
    if (nImplicitCalls < nSampleImpl)
    {
        cout << "could not get enough number of samples for implicit solver" << endl;
        cout << "expected number of samples = " << nSampleImpl << endl;
        cout << "real number of samples = " << nImplicitCalls << endl;
    }
    else
    {
        cout << "time spent for implicit solver for " << nSampleImpl << " samples = " << fixed << setprecision(2) << (eTimeImpl /= nSampleImpl) << wImpl.unit << endl;
    }
    
}