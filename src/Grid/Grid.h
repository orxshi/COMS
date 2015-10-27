/* 
 * File:   Grid.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 5:41 AM
 */

#ifndef GRID_H
#define	GRID_H

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <functional>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include "../Vector/Vector.h"
#include "../Face/Face.h"
#include "BinaryTree.h"
#include "../ADT/ADT.h"
#include "../LinearAlgebra/LinearAlgebra.h"
#include "../Point/Point.h"
#include "../Solid/Solid.h"

using std::ifstream;
using std::ofstream;
using std::vector;
using std::string;
using std::ref;
using std::cout;
using std::endl;
using std::exit;
using std::move;
using std::to_string;
using std::min;
using std::max;
using std::setprecision;
using std::fixed;
using std::function;
using std::flush;
using std::left;
using std::setw;
using std::stringstream;
using std::scientific;

struct Hole
{
    array <double, N_DIM> min;
    array <double, N_DIM> max;
};

struct Grid
{
    // Fields
    int n_bou_elm;
    int n_in_elm;
    int totalNElms;
    int id;
    int phys_count;
    int nHoles;    
    double wallDistance;    
    ifstream in; // for reading grid input
    string meshFile;
    string outputDir;
    string mainDir;
    string logDir;
    vector<int> phys;
    vector<int> bc;
    vector <Point> pt;
    vector <Face> face;
    vector <Cell> cell;
    vector <btree> bt;
    vector <Hole> holes;    
    array<string,3> bcVerbose;    
    
    struct CellADT: ADT
    {
        bool compareFunction (const Node *node, const ADTPoint& targetPoint) override;
        bool doCubesOverlap (const Node *node, const ADTPoint& targetPoint) override;
        void build (const Grid& gr);
    } cellADT;

    // Constructor
    Grid (string mainDir, int id);
    //Grid (const Grid&);
    Grid (Grid&& other);
    Grid& operator=(Grid&& other);
    
    // Destructor
    ~Grid ();

    // Private Methods
    void read_ptSize ();
    void read_pt ();
    void read_elmSize ();
    void read_elm ();
    void set_faceVertices (Face& face, const Cell& elm, const int index);
    void set_connectivity ();
    void set_elmVolumes ();
    void set_elmCentroids ();
    void updateVars();
    double setExpRes();
    void printInput();
    void printMeshInfo();

    // Public Methods
    void readInput();    
    void read_grid();
    void set_grid();
    void apply_BCs();
    void set_BCs();
    //void createOutputDir(const string mainDir);
    void identifyIBlank(Grid& gr);
    void leastSquaresCoeffs();
    void leastSquaresGrad();
    void interpolate();
    void trimWhoHasFringeNeighbor();
    void trimWhoHasTrimNeighbor (int threshold);
    void trimToUntrim (int crt);
    void nontrimToTrim (int crt);
    void fieldToFringe (int crt);
    void fringeToField (int crt);
    //void blankWithPhys (int phys);
    void getMovingFaceVelocity();
    void rotateMesh (double delAlpha, CVector& centerAirfoil);    
    void outAllTecplot();
    void outAllUnsteady(string title);
    void log();
    void setWallDistance (int phys);
    void outAllVTK (int time);    
};

#endif	/* GRID_H */

