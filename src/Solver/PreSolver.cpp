#include "../Solver/Solver.h"

// init

/*void Solver::preSolver(Grid& gr)
{
    gr.init();
    cout << ">>> Setting BCs..." << flush;
    gr.set_BCs();
    cout << ">>> done!" << endl;
    cout << ">>> Applying BCs..." << flush;
    gr.apply_BCs();
    cout << ">>> done!" << endl;
    
    for (Cell& cll: gr.cell)
    {
        cll.old_cons = cll.cons;
    }
    
    if (sOrder == tsOrder_t::SECOND)
    {
        gr.leastSquaresCoeffs();
    }
}*/
