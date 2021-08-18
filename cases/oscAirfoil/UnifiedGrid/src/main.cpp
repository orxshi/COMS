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
    
    Watch watchSteady;
    Watch watchOscAirfoil;
    
    //PetscMPIInt rank, n_procs;
    //MPI_Comm world = PETSC_COMM_WORLD;
    //MPI_Comm_rank (world, &rank);
    //MPI_Comm_size (world, &n_procs);
    
    string mainDir = createOutputDir();
    
    // background grid
    Grid bg (mainDir, 0);
    bg.read_grid();
    bg.set_grid();    
    
    // airfoil grid
    Grid ag (mainDir, 1);
    ag.read_grid();
    ag.set_grid();    
    
    // initialize grids
    OscInit oscInit;
    oscInit.read();
    oscInit.init (ag);
    oscInit.init (bg);
    
    // push grids to vector
    vector<Grid> grs;    
    grs.push_back(move(bg));
    grs.push_back(move(ag));
    
    // set wall distances
    grs[0].setWallDistance(2);
    grs[1].setWallDistance(3);
    
    //grs[0].cellADT.build (grs[0]);
    //grs[1].cellADT.build (grs[1]);
    //Iblank iBlank;
    //iBlank.identify (grs[0], grs[1]);
    //iBlank.identify (grs[1], grs[0]);

    
    //grs[0].outAllVTK (0);
    //grs[1].outAllVTK (0);
    
    Grid finalGrid (mainDir, 3);

    
    //AFT::aft (grs, finalGrid);    
    //finalGrid.outAllVTK (0);

    //return 0;

    
    //finalGrid.readInput();
    //finalGrid.leastSquaresCoeffs();    
    //finalGrid.cellADT.build (finalGrid);
    //oscInit.init (finalGrid);
    
    //Solver solSteady (finalGrid, "SOLVER-STEADY");
    //solSteady.read ("Solver/solSteady.dat");

    // solve steady state
    //SMAirfoil sma (solSteady.dt);
    OscAirfoil oa (1.); // 1 is time step
    //sma.read ("MovingGrid/smAirfoil.dat");
    oa.read ("MovingGrid/oscAirfoil.dat");

    //Coeffs coeffs (finalGrid, oscInit.rhoInf, oscInit.pInf, oscInit.Mach, oa.MachAirfoil);
    
    //sma.getAllFaceVelocities (finalGrid);
    //watchSteady.start();
    //(solSteady.implicit) ? solSteady.impl(finalGrid) : solSteady.expl(finalGrid);
    //solSteady.petsc.finalize();
    //watchSteady.stop();
    
    //finalGrid.outAllVTK (0);
    //exit(-2);
        
    // solve osc airfoil
    //Grid oldGrid (mainDir, 4);
    //Grid oldGrid = move(finalGrid);
    int countr = 0;
    watchOscAirfoil.start();
    //for (double time=0.; time<100.; time+=1.0) // 1 is dt
    for (double time=0.; time<1.; time+=1.0) // 1 is dt
    {
//centerAirfoil[0] 0.25
//centerAirfoil[1] 0.
//centerAirfoil[2] 0.
double kc = 0.0814;
double MachAirfoil = -0.75;
double alphaMean = 0.016;
double alphaMax = 2.51;
        double tmp = 2. * kc * fabs(MachAirfoil);

        double alpha    = alphaMean + alphaMax * sin( tmp * time ); // degree
        double delAlpha = -tmp * alphaMax * cos( tmp * time ); // degree
        cout << "time = " << time << endl;
        cout << "alpha = " << alpha << endl;

        grs[0].outAllVTK (countr);
        grs[1].outAllVTK (countr);
        
        grs[0].cellADT.build (grs[0]);
        grs[1].cellADT.build (grs[1]);   
        Iblank iBlank;
        iBlank.identify (grs[0], grs[1]);
        iBlank.identify (grs[1], grs[0]);

        Grid finalGrid (mainDir, 3);

        AFT::aft (grs, finalGrid, countr);
        //finalGrid.cellADT.build (finalGrid);        
        //finalGrid.readInput();
        //finalGrid.leastSquaresCoeffs();
        
        if (time == 0.)
        {            
            /*oa.delAlpha = 0.;
            oscInit.init (finalGrid);
            oa.interFromOldTS (finalGrid, oldGrid);*/
            //finalGrid = move(oldGrid);
        }
        else
        {   
            /*oscInit.init (finalGrid);
            oa.interFromOldTS (finalGrid, oldGrid);*/
            //finalGrid.set_BCs();
            //finalGrid.apply_BCs();
        }
        
        //Solver solOscAirfoil (finalGrid, "SOLVER-OSC-AIRFOIL");
        //solOscAirfoil.read ("Solver/solOscAirfoil.dat");
        //solOscAirfoil.time = time;
        
        oa.setAngles (time);
        /*oa.getAllFaceVelocities (finalGrid);
        (solOscAirfoil.implicit) ? solOscAirfoil.impl(finalGrid) : solOscAirfoil.expl(finalGrid);
        coeffs.getCoeffs (finalGrid);
        outLiftCoef (coeffs, oa.alpha, solOscAirfoil.time);
        coeffs.outPresCoef (countr);*/
        finalGrid.outAllVTK (countr);
        oa.moveGrid (grs[1]);
        
        //oldGrid = move(finalGrid);
        
        ++countr;
    }
    watchOscAirfoil.stop();
    
    /*if (rank == MASTER_RANK)
    {
        //gr.outAllTecplot();
        finalGrid.outAllVTK (0);
        //coeffs.out.close();
        log (mainDir, watchSteady.elapsedTime, "elapsedTimeSteady", watchSteady.unit);
        //log (mainDir, watchOscAirfoil.elapsedTime, "elapsedTimeOscAirfoil", watchOscAirfoil.unit);
        solSteady.log (finalGrid.logDir);
        //solOscAirfoil.log (gr.logDir);
        sma.log (finalGrid.logDir);
        //oa.log (gr.logDir);
    }*/
    
    //PetscFinalize();

    return 0;
}
