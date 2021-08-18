#include "Solver.h"

inline void invertMat5(const Matrixd<5,5>& A, const Vector<N_VAR>& f, Vector<N_VAR>& x)
{
    // perform x = inv (A)*f

    double b11,b21,b22,b31,b32,b33,b41,b42,b43,b44,b51,b52,b53,b54,b55;
    double u12,u13,u14,u15,u23,u24,u25,u34,u35,u45;
    double d1,d2,d3,d4,d5;

    // decompose A into L and U
    b11 = 1 / A(0,0);
    u12 = A(0,1) * b11;
    u13 = A(0,2) * b11;
    u14 = A(0,3) * b11;
    u15 = A(0,4) * b11;
    b21 = A(1,0);
    b22 = 1 / (A(1,1) - b21 * u12);
    u23 = (A(1,2) - b21 * u13) * b22;
    u24 = (A(1,3) - b21 * u14) * b22;
    u25 = (A(1,4) - b21 * u15) * b22;
    b31 = A(2,0);
    b32 = A(2,1) - b31 * u12;
    b33 = 1 / (A(2,2) - b31 * u13 - b32 * u23);
    u34 = (A(2,3) - b31 * u14 - b32 * u24) * b33;
    u35 = (A(2,4) - b31 * u15 - b32 * u25) * b33;
    b41 = A(3,0);
    b42 = A(3,1) - b41 * u12;
    b43 = A(3,2) - b41 * u13 - b42 * u23;
    b44 = 1 / (A(3,3) - b41 * u14 - b42 * u24 - b43 * u34);
    u45 = (A(3,4) - b41 * u15 - b42 * u25 - b43 * u35) * b44;
    b51 = A(4,0);
    b52 = A(4,1) - b51 * u12;
    b53 = A(4,2) - b51 * u13 - b52 * u23;
    b54 = A(4,3) - b51 * u14 - b52 * u24 - b53 * u34;
    b55 = 1 / (A(4,4) - b51 * u15 - b52 * u25 - b53 * u35 - b54 * u45);
    //
    d1 = f[0] * b11;
    d2 = (f[1] - b21 * d1) * b22;
    d3 = (f[2] - b31 * d1 -b32 * d2) *b33;
    d4 = (f[3] - b41 * d1 -b42 * d2 - b43 * d3) * b44;
    d5 = (f[4] - b51 * d1 -b52 * d2 - b53 * d3 - b54 * d4) * b55;
    //
    x[4] = d5;
    x[3] = d4 - u45 * d5;
    x[2] = d3 - u34 * x[3] - u35 * d5;
    x[1] = d2 - u23 * x[2] - u24 * x[3] - u25 * d5;
    x[0] = d1 - u12 * x[1] - u13 * x[2] - u14 * x[3] - u15 * d5;
}

/*inline Vector<MAT5_SIZE> mat5Vec5Mul (const Matrix5& M, const Vector<MAT5_SIZE>& V)
{
    Vector<MAT5_SIZE> T;
    T.fill(0.);
    
    for (int r=0; r<MAT5_SIZE; ++r)
    {
        for (int c=0; c<MAT5_SIZE; ++c)
        {
            T[r] += M[r][c] * V[c];
        }
    }
    
    return T;
}*/

inline void common2 (Cell& e, const vector<Face>& face, vector<Cell>& cell)
{
    Vector<N_VAR> sum;
    sum[0] = 0.;
    sum[1] = 0.;
    sum[2] = 0.;
    sum[3] = 0.;
    sum[4] = 0.;
    
    for (int i=0; i<e.nei.size(); ++i)
    {
        const Face& f = face[e.face[i]];
        Cell& LC = cell[ f.nei[0] ];
        Cell& RC = cell[ f.nei[1] ];
        
        if (&e == &LC)
        {
            //sum -= mat5Vec5Mul(f.M[1], cell[e.nei[i]].dQ);
            //sum -= f.M[1] % cell[e.nei[i]].dQ;
        }
        else
        {
            //sum += mat5Vec5Mul(f.M[0], cell[e.nei[i]].dQ);
            //sum += f.M[0] % cell[e.nei[i]].dQ;
        }
    }
    
    sum = sum + e.R;
    //invertMat5(e.D, sum, e.dQ);
}

inline void common (Cell& e, const vector<Face>& face, vector<Cell>& cell, const vector <Matrixd<N_VAR,N_VAR>>& M0, const vector <Matrixd<N_VAR,N_VAR>>& M1)
{
    invertMat5(e.D, e.R, e.dQ); // dQ = inv(D) * res

    //std::cout << "R" << std::endl;
    //std::cout << e.R[0] << " " << e.R[1] << " " << e.R[2] << " " << e.R[3] << " " << e.R[4] << " " << std::endl;

    //std::cout << "D" << std::endl;
    //std::cout << e.D(0,0) << " " << e.D(0,1) << " " << e.D(0,2) << " " << e.D(0,3) << " " << e.D(0,4) << " " << std::endl;
    //std::cout << e.D(1,0) << " " << e.D(1,1) << " " << e.D(1,2) << " " << e.D(1,3) << " " << e.D(1,4) << " " << std::endl;
    //std::cout << e.D(2,0) << " " << e.D(2,1) << " " << e.D(2,2) << " " << e.D(2,3) << " " << e.D(2,4) << " " << std::endl;
    //std::cout << e.D(3,0) << " " << e.D(3,1) << " " << e.D(3,2) << " " << e.D(3,3) << " " << e.D(3,4) << " " << std::endl;
    //std::cout << e.D(4,0) << " " << e.D(4,1) << " " << e.D(4,2) << " " << e.D(4,3) << " " << e.D(4,4) << " " << std::endl;

    //std::cout << "dQ" << std::endl;
    //std::cout << e.dQ[0] << " " << e.dQ[1] << " " << e.dQ[2] << " " << e.dQ[3] << " " << e.dQ[4] << " " << std::endl;

    //assert(false);
    
    Vector<N_VAR> ddQ = e.dQ - e.old_dQ;

    // update neighbour residuals
    for (const int i: e.face)
    {
        const Face& f = face[i];
        Cell& LC = cell[ f.nei[0] ];
        Cell& RC = cell[ f.nei[1] ];
        
        if (&e == &LC)
        {
            //RC.R += mat5Vec5Mul(f.M[0], ddQ);
            //RC.R += f.M[0] % ddQ;
            RC.R += M0[i] % ddQ;
        }
        else
        {
            //LC.R -= mat5Vec5Mul(f.M[1], ddQ);
            //LC.R -= f.M[1] % ddQ;
            LC.R -= M1[i] % ddQ;
        }
    }
}

void Solver::diff_to_cons_prim(Grid& g)
{
    for (int ic=g.n_bou_elm; ic<g.cell.size(); ++ic)
    {
        Cell& e = g.cell[ic];
        
        if (e.iBlank == iBlank_t::FIELD)
        {
            e.cons += e.dQ;
            //if (e.dQ[0] > 1)
            //{
                //assert(false);
            //}
            e.cons_to_prim();
        }
    }
}

void Solver::set_residual(Grid& g)
{
    Vector<N_VAR> res;
    //res.fill(0.);
    res.fill(BIG_NEG_NUM);
    
    for (int ic=g.n_bou_elm; ic<g.cell.size(); ++ic)
    {
        Cell& e = g.cell[ic];
        
        e.sigma = 0.;
        e.R.fill(0.);
        //eq5 (e.D, 0.);
        e.D = 0.;
    }    
    //if (sOrder == 2) { g.leastSquaresGrad(); }
    //roeflx (g);
    
    for (int ic=g.n_bou_elm; ic<g.cell.size(); ++ic)
    {
        Cell& e = g.cell[ic];
        
        if (e.iBlank == iBlank_t::FIELD)
        {
            for (auto i=0; i<e.dQ.size(); ++i)
            {
                if ( std::isnan(e.dQ[i]) ) { cout << "e.dQ[i] is NAN in Solver::set_residual(...)" << endl; exit(-2); }
                res[i] = max(fabs(e.R[i]),res[i]);
                //res[i] += pow( e.dQ[i],2 );
                
                /*if (e.dQ[i] != e.R[i])
                {
                    cout << "dsfgdfgdfg" << endl;
                    cout << "dQ[i] = " << e.dQ[i] << endl;
                    cout << "R[i] = " << e.R[i] << endl;
                }*/
            }
        }
    }
    
    aveRes = BIG_NEG_NUM;
    for (auto i=0; i<N_VAR; ++i)
    {
        //res[i] /= g.n_in_elm;
        //res[i] = sqrt(res[i]);
        aveRes = max(res[i],aveRes);
    }
    
    //aveRes = std::accumulate (res.begin(), res.end(),0.);
    //aveRes /= N_VAR;
    
    if ( std::isnan(aveRes) ) { cout << "aveRes is NAN in Solver::set_residual(...)" << endl; exit(-2); }
    if ( std::isinf(aveRes) ) { cout << "aveRes is INF in Solver::set_residual(...)" << endl; exit(-2); }
}

void Solver::gauss_seidel (Grid& g)
{
    for (int ic=g.n_bou_elm; ic<g.cell.size(); ++ic)
    {
        g.cell[ic].old_dQ.fill(0.);
    }
    
    for (int i=0; i<nGaussIter; ++i)
    {
        // forward
        for (int ic=g.n_bou_elm; ic<g.cell.size()-1; ++ic)
        {          
            Cell& e = g.cell[ic];
        
            if (e.iBlank == iBlank_t::FIELD)
            {
                common (e, g.face, g.cell, M0, M1);
                e.old_dQ = e.dQ;
            }
        }

        // backward
        for (int ic=g.cell.size()-1; ic>=g.n_bou_elm; --ic)
        {
            Cell& e = g.cell[ic];
        
            if (e.iBlank == iBlank_t::FIELD)
            {
                common (e, g.face, g.cell, M0, M1);
                e.old_dQ = e.dQ;
            }
        }
    }
}
