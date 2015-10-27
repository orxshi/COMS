#include "Method.h"

/*void presetMulti (vector<Grid>& gr)
{
    string mainDir = createOutputDir();
    
    unsigned int size = gr.size();
    
    for (int i=0; i<size; ++i)
    {
        gr[i].id = i;
        gr[i].read_grid (i);
        gr[i].set_grid();
        gr[i].createOutputDir (mainDir);
        gr[i].preSolver();
    }
    
    for (int i=0; i<size; ++i)
    {
        cout << ">>> Building grid[" << i << "].cellADT" << endl;

        gr[i].setWallDistance();
        gr[i].cellADT.build (gr[i]);

        cout << ">>> Built grid[" << i << "].cellADT" << endl;
    }

    cout << ">>> Identifying (non)active cells" << endl;

    for (int i=0; i<size; ++i)
    {
        for (int j=0; j<size; ++j)
        {
            if (j != i)
            {
                gr[i].identifyIBlank (gr[j]);
            }
        }
    }

    cout << ">>> Finished identifying (non)active cells" << endl;
}*/

/*void presetSingle (Grid& gr)
{
    string mainDir = createOutputDir();
    
    gr.id = 0;
    gr.readInput();
    gr.read_grid();
    gr.set_grid();
    gr.createOutputDir (mainDir);
    gr.preSolver();
    
    for (Cell& cll: gr.cell)
    {
        cll.iBlank = iBlank_t::FIELD;
    }
}*/

/*void postset (vector<Grid>& gr)
{
    for (Grid& g: gr)
    {        
        g.outAllTecplot();
        g.logFile();
    }
}*/

/*void SSG (Grid& gr, Solver& sol)
{
    // (S)olve (S)ingle (G)rid
    
    (sol.implicit) ? sol.impl(gr) : sol.expl(gr);
}*/

/*void SOGTM (vector<Grid>& gr, double tol)
{
    // (S)olve (O)verset (G)rids with (T)raditional (M)ethod
    
    unsigned int convCounter;
    unsigned int size = gr.size();
    bool verbose = true;
    double res;
    
    do
    {
        convCounter = 0;

        for (unsigned int i=0; i<size; ++i)
        {
            gr[i].interpolate();
            gr[i].apply_BCs();
            res = gr[i].setExpRes();
            cout << "Grid_" << i << "_res = " << res << endl;
            if (res < tol) { ++convCounter; }
            (gr[i].implicit) ? gr[i].impl(verbose) : gr[i].expl(verbose);
        }
        
        cout << "convCounter = " << convCounter << endl;
    }
    while  ( convCounter != size );
}*/

/*void SOGNM (vector<Grid>& gr)
{
    // (S)olve (O)verset (G)rids with (N)ew (M)ethod
    
    Grid finalGrid;
    
    AFT::aft (gr, finalGrid);
    
    finalGrid.preSolver();
    SSG (finalGrid);
    finalGrid.outAllTecplot();
    finalGrid.logFile();
}*/

/*void customUnsteady (vector<Grid>& gr)
{
    unsigned int size = gr.size();
    
    for (unsigned int i=0; i<size; ++i)
    {
        gr[i].getMMVariables();
        gr[i].getMovingFaceVelocity();
        gr[i].apply_BCs();
    
        gr[i].alpha = 0.016;
        gr[i].delAlpha = 0.;
        gr[i].globalStep = -1;
    }
    
    SSG (gr[0]);
    
    for (unsigned int i=0; i<size; ++i)
    {
        gr[i].steady = false;
        gr[i].useCFL = false;
        //gr[i].cfl = 10;
        gr[i].kc = 0.0814;
        gr[i].alphaMax = 2.51;
        gr[i].alphaMean = 0.016;
    }
    
    double dt = gr[0].dt;
    gr[0].globalStep = 0;
    
    for (double time=0.; time<52; time += dt)
    {
        //double startTime = getWallTime();
        
        for (unsigned int i=0; i<size; ++i)
        {
            gr[i].time = time;
        }
        
        for (unsigned int i=0; i<size; ++i)
        {
            gr[i].getMMVariables();
            gr[i].getMovingFaceVelocity();
            gr[i].apply_BCs();
        }
        
        SSG (gr[0]);
        
        for (unsigned int i=0; i<size; ++i)
        {
            gr[i].coeffs();
            gr[i].rotateMesh();
        }
        gr[0].outAllUnsteady();
    }
}*/

/*void chooseMethod(method_t& method)
{
    int nGrids;
    int tmpi;
    string tmps;
    ifstream is;
    double oversetTol;
    
    is.open ("../input.dat");
    is >> tmps; is >> tmpi;
    
    switch (tmpi)
    {
        case 0:
            method = method_t::SINGLE;
            cout << ">>> Method: SINGLE" << endl;
            break;
        case 1:
            method = method_t::MULTI;
            cout << ">>> Method: MULTI" << endl;
            break;
        case 2:
            method = method_t::NEW;
            cout << ">>> Method: NEW" << endl;
            break;
        default:
            method = method_t::UNDEFINED;
            cout << ">>> Method: UNDEFINED" << endl;
            break;
    }
    
    if (method == method_t::MULTI || method == method_t::NEW)
    {
        is >> tmps; is >> nGrids;
        
        if ( method == method_t::MULTI )
        {
            is >> tmps; is >> oversetTol;
        }
    }
    else if (method == method_t::SINGLE)
    {
        nGrids = 1;
    }
    else
    {
        cout << "unknown method" << endl;
        exit(-2);
    }
    
    cout << ">>> nGrids: " << nGrids << endl;
    cout << string(SPLITTER_LEN, '-') << endl;
    
    vector<Grid> gr (nGrids);
    preset (gr);
}*/

/*void chooseMethod()
{
    double startTime, endTime;
    
    chooseMethod (method);
    
    switch (method)
    {
        case method_t::SINGLE:
            startTime = getWallTime();
            //SSG (gr[0]);
            customUnsteady (gr);
            endTime = getWallTime();
            getElapsedTime (startTime, endTime, gr[0].elapsedTime, gr[0].elapsedTimeUnit);
            postset (gr);
            break;
        case method_t::MULTI:
            startTime = getWallTime();
            //SOGTM (gr, oversetTol);
            customUnsteady (gr);
            endTime = getWallTime();
            for (Grid& g: gr)
            {
                getElapsedTime (startTime, endTime, g.elapsedTime, g.elapsedTimeUnit);
            }
            postset (gr);
            break;
        case method_t::NEW:
            startTime = getWallTime();
            SOGNM (gr);
            endTime = getWallTime();
            getElapsedTime (startTime, endTime, gr[0].elapsedTime, gr[0].elapsedTimeUnit);
            break;
        default:
            cout << "unknown method" << endl;
            exit(-2);
            break;
    }
}*/
