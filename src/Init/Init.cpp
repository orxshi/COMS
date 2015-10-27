#include "../Grid/Grid.h"

/*void Grid::init ()
{
    //double rhoInf, pInf, Mach, aoaDeg;
    //string tmps;    
    //ifstream in;
    //in.open("./Init/init.dat");
    
    /*in >> tmps; in >> rhoInf;
    in >> tmps; in >> pInf;
    in >> tmps; in >> Mach;
    in >> tmps; in >> aoaDeg;*/
    
    /*cout << ">>> Ref density: " << rhoInf << endl;
    cout << ">>> Ref pressure: " << pInf << endl;
    cout << ">>> Mach: " << Mach << endl;
    cout << ">>> aoa (deg): " << aoaDeg << endl;
    cout << ">>> done!" << endl;
    //cout << string(SPLITTER_LEN, '-') << endl;*/
    
    /*double aoaRad = aoa * DEG_TO_RAD;
    
    double p[5];
    
    p[0] = rhoInf;
    p[1] = MachAir * cos(aoaRad);
    p[2] = MachAir * sin(aoaRad);
    p[3] = 0.;
    p[4] = pInf;

    for (Cell& cll: cell)
    {
        for (int i=0; i<N_VAR; ++i)
        {
            cll.prim[i] = p[i];
        }
        
        cll.prim_to_cons();
    }
}*/