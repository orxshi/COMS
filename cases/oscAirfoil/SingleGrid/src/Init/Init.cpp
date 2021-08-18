#include "Init.h"

void OscInit::init (Grid& gr)
{    
    double aoaRad = aoa * DEG_TO_RAD;
    
    double p[5];
    
    p[0] = rhoInf;
    p[1] = Mach * cos(aoaRad);
    p[2] = Mach * sin(aoaRad);
    p[3] = 0.;
    p[4] = pInf;

    //std::cout << "rhoInf: " << rhoInf << std::endl;

    for (Cell& cll: gr.cell)
    {
        for (int i=0; i<N_VAR; ++i)
        {
            cll.prim[i] = p[i];
        }
        
        cll.prim_to_cons();
        
        cll.old_cons = cll.cons;
        cll.oldold_cons = cll.cons;
        
        cll.iBlank = iBlank_t::FIELD;
    }
    
    gr.set_BCs();
    gr.apply_BCs();
}

void OscInit::init_sod (Grid& gr)
{    
    double pL[5], pR[5];
    
    pL[0] = 1.;
    pL[1] = 0.;
    pL[2] = 0.;
    pL[3] = 0.;
    pL[4] = 100;

    pR[0] = 0.125;
    pR[1] = 0.;
    pR[2] = 0.;
    pR[3] = 0.;
    pR[4] = 10;

    for (Cell& cll: gr.cell)
    {
        if (cll.cnt[0] < 0.)
        {
            for (int i=0; i<N_VAR; ++i)
            {
                cll.prim[i] = pL[i];
            }
        }
        else
        {
            for (int i=0; i<N_VAR; ++i)
            {
                cll.prim[i] = pR[i];
            }
        }
        
        cll.prim_to_cons();
        
        cll.old_cons = cll.cons;
        cll.oldold_cons = cll.cons;
        
        cll.iBlank = iBlank_t::FIELD;
    }
    
    gr.set_BCs();
    gr.apply_BCs();
}
