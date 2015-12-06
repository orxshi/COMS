#include "../Solver/Solver.h"

void Solver::interflux (Grid& gr)
{
    double dTao;
    
    for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
    {
        Cell& e = gr.cell[ic];
        
        dTao = cfl * e.vol / e.sigma;
        
        switch (tOrder)
        {
            case 1:
                if (steady)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        e.D(i,i) += e.vol / dTao;
                    }
                }
                else
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        e.D(i,i) += e.vol * (1./dTao + 1./dt);
                        e.R[i] -= e.vol * (e.cons[i] - e.old_cons[i]) / dt;
                    }
                }
                break;
            case 2:
                if (steady)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        e.D(i,i) += 1.5 * e.vol / dTao;
                    }
                }
                else
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        e.D(i,i) += 1.5 * e.vol * (1./dTao + 1./dt);
                    }
                    
                    if (nTimeStep > 1)
                    {
                        for (int i=0; i<N_VAR; ++i)
                        {
                            e.R[i] -= 0.5 * e.vol * (3.*e.cons[i] - 4.*e.old_cons[i] + e.oldold_cons[i]) / dt;
                        }
                    }
                }
                break;
            default:
                cout << "undefined tOrder in Solver::interflux(...)" << endl;
                exit(-2);
                break;
        }
    }
}
