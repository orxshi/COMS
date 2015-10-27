#include <petscsys.h>
#include "../../../../src/Grid/Grid.h"
#include "../../../../src/Output/Output.h"
#include "../../../../src/Solid/Solid.h"
#include "../../../../src/Coeffs/Coeffs.h"
#include "../../../../src/Solver/Solver.h"
#include "../../../../src/AFT/AFT.h"
#include "Init/Init.h"
#include "MovingGrid/OscAirfoil.h"
#include "MovingGrid/StraightMovingAirfoil.h"
#include "Output/Output.h"

int main(int argc, char** argv)
{
    PetscInitialize (&argc, &argv, NULL, NULL);
    
    Watch watchSteady;
    Watch watchOscAirfoil;
    
    PetscMPIInt rank, n_procs;
    MPI_Comm world = PETSC_COMM_WORLD;
    MPI_Comm_rank (world, &rank);
    MPI_Comm_size (world, &n_procs);
    
    string mainDir = createOutputDir();
    
    // background grid
    Grid bg (mainDir, 1);
    bg.read_grid();
    bg.set_grid();    
    
    // airfoil grid
    Grid ag (mainDir, 0);
    ag.read_grid();
    ag.set_grid();    
    
    // initialize grids
    OscInit oscInit;
    oscInit.read();
    oscInit.init (ag);
    oscInit.init (bg);
    
    // push grids to vector
    vector<Grid> grs;    
    grs.push_back(move(ag));
    grs.push_back(move(bg));
    
    // set wall distances
    grs[0].setWallDistance(3);
    grs[1].setWallDistance(2);
    
    grs[0].cellADT.build (grs[0]);
    grs[1].cellADT.build (grs[1]);        
    grs[0].identifyIBlank (grs[1]);
    grs[1].identifyIBlank (grs[0]);
    
    array<Solver,2> solver = { Solver(grs[0], "SOLVER-STEADY-AG"), Solver(grs[1], "SOLVER-STEADY-BG") };    
    solver[0].read ("Solver/solSteady.dat");
    solver[1].read ("Solver/solSteady.dat");
    
    SMAirfoil sma (solver[0].dt);
    sma.read ("MovingGrid/smAirfoil.dat");
    
    double err = BIG_POS_NUM;
    //while (err > solver[0].tol)
    while (false)
    {
        for (int g=0; g<grs.size(); ++g)
        {
            sma.getAllFaceVelocities (grs[g]);
            
            (solver[g].implicit) ? solver[g].impl(grs[g]) : solver[g].expl(grs[g]);
            
            grs[1-g].interpolate();
            sma.getAllFaceVelocities (grs[1-g]);
            
            err = solver[1-g].setExpRes (grs[1-g]);
        }
        
        cout << "err = " << err << endl;
    }
    
    solver[0].petsc.finalize();
    solver[1].petsc.finalize();
    
    grs[0].outAllVTK (0);
    grs[1].outAllVTK (0);
    
    PetscFinalize();

    return 0;
}