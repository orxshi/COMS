#include "Solver.h"

void Solver::getRes (Grid& gr, Limiter& limiter)
{
    Vector<N_VAR> res;
    //res.fill(0.);    
    res.fill(BIG_NEG_NUM); // for max criterion
    
    // calculate cll.R with updated values. updateVars should be called before this
    
    if (sOrder == 2)
    {
        //gr.leastSquaresGrad();
        gradient.leastSquaresGrad (gr);
    }
    
    roe.roeflx (gr, limiter, M0, M1, gradient);
    
    if (tOrder == 1)
    {
        if (steady)
        {
            for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
            {
                Cell& cll = gr.cell[ic];
            
                if (cll.iBlank == iBlank_t::FIELD)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        double RHS = cll.R[i];
                        res[i] = max( fabs(RHS), res[i] );
                    }
                }
            }
        }
        else
        {
            for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
            {
                Cell& cll = gr.cell[ic];
            
                if (cll.iBlank == iBlank_t::FIELD)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        double LHS = (cll.cons[i] - cll.old_cons[i]) * cll.vol / dt;
                        double RHS = cll.R[i];

                        //res[i] += pow(LHS - RHS, 2.);
                        res[i] = max( fabs(LHS - RHS), res[i] );
                    }
                }
            }
        }
    }
    else
    {
        cout << "time derivative is not 1 or 2 in Solver::setExpRes(...)" << endl;
        exit(-2);
    }
    
    aveRes = BIG_NEG_NUM;
    for (int i=0; i<N_VAR; ++i)
    {
        aveRes = max(res[i],aveRes);
    }
    
    if ( isnan(aveRes) ) { cout << "averageRes is NAN in setExpRes()" << endl; exit(-2); }
    if ( isinf(aveRes) ) { cout << "averageRes is INF in setExpRes()" << endl; exit(-2); }
}