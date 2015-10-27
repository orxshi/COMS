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