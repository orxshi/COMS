//#include <petscsys.h>
#include "../../../../src/Grid/Grid.h"
#include "../../../../src/IBlank/IBlank.h"
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
    //PetscInitialize (&argc, &argv, NULL, NULL);
    //MPI_Init(NULL, NULL);
    
    Watch watchSteady;
    Watch watchOscAirfoil;
    
    //PetscMPIInt rank, n_procs;
    //MPI_Comm world = PETSC_COMM_WORLD;
    //MPI_Comm_rank (world, &rank);
    //MPI_Comm_size (world, &n_procs);
    
    string mainDir = createOutputDir();
    vector<Grid> grs = {Grid(mainDir, 1), Grid(mainDir, 0)};    
    
    // background grid
    //Grid bg (mainDir, 1);
    grs[0].read_grid();
    grs[0].set_grid();    
    
    // airfoil grid
    //Grid ag (mainDir, 0);
    grs[1].read_grid();
    grs[1].set_grid();    
    
    // initialize grids
    OscInit oscInit;
    oscInit.read();
    oscInit.init (grs[1]);
    oscInit.init (grs[0]);
    
    // push grids to vector
    //vector<Grid> grs;    
    //grs.push_back(move(ag));
    //grs.push_back(move(bg));
    
    // set wall distances
    grs[0].setWallDistance(3);
    grs[1].setWallDistance(2);
    
    // build cell trees
    grs[0].cellADT.build (grs[0]);
    grs[1].cellADT.build (grs[1]);        
    
    // hole cutting
    Iblank iblank;
    iblank.identify (grs[0], grs[1]);
    iblank.identify (grs[1], grs[0]);
    
    grs[0].outAllVTK (0);
    grs[1].outAllVTK (0);

    //// solvers
    //array<Solver,2> solver = { Solver(grs[0], "SOLVER-STEADY-AG"), Solver(grs[1], "SOLVER-STEADY-BG") };    
    //solver[0].read ("Solver/solSteady.dat");
    //solver[1].read ("Solver/solSteady.dat");
    //
    //// moving airfoil object
    //SMAirfoil sma (solver[0].dt);
    //sma.read ("MovingGrid/smAirfoil.dat");
    //
    //int nActiveElms = 0;
    //for (int g=0; g<grs.size(); ++g)
    //{
    //    for (const Cell& cll: grs[g].cell)
    //    {
    //        if (cll.iBlank == iBlank_t::FIELD)
    //        {
    //            ++nActiveElms;
    //        }
    //    }
    //}
    //log (mainDir, nActiveElms, "nActiveElms", "");
    //cout << "nActiveElms = " << nActiveElms << endl;
    //
    //watchSteady.start();
    //double err0 = BIG_POS_NUM;
    //double err1 = BIG_POS_NUM;
    //solver[0].maxTimeStep = 1;
    //solver[1].maxTimeStep = 1;
    //double oversetTol = 1e-12;
    //while (err0 > oversetTol || err1 > oversetTol)    
    //{
    //    for (int g=0; g<grs.size(); ++g)
    //    {
    //        sma.getAllFaceVelocities (grs[g]);
    //        
    //        (solver[g].implicit) ? solver[g].impl(grs[g]) : solver[g].expl(grs[g]);
    //        
    //        //grs[1-g].interpolate();
    //        iblank.interpolate (grs[1-g], solver[g].gradient);

    //        sma.getAllFaceVelocities (grs[1-g]);
    //        
    //        err0 = solver[1-g].setExpRes (grs[1-g]);
    //        err1 = solver[g].setExpRes (grs[g]);
    //    }
    //    
    //    cout << "err0 = " << err0 << endl;
    //    cout << "err1 = " << err1 << endl;
    //}
    //watchSteady.stop();
    //log (mainDir, watchSteady.elapsedTime, "elapsedTimeSteady", watchSteady.unit);
    //
    ////solver[0].petsc.finalize();
    ////solver[1].petsc.finalize();
    //
    //grs[0].outAllVTK (0);
    //grs[1].outAllVTK (0);
    //
    ////PetscFinalize();
    //MPI_Finalize();

    return 0;
}
