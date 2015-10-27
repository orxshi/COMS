#include "Roe.h"

Matrix5 jacob(const Vector<N_VAR>& q, const CVector& n, double vbn)
{
    Matrix5 A;
    
    double gs = GAMMA - 1.;

    double nx = n[0];
    double ny = n[1];
    double nz = n[2];
    
    double kq = 0.5 * ( pow(q[1],2) + pow(q[2],2) + pow(q[3],2) );
    double qn = q[1]*nx + q[2]*ny + q[3]*nz;

    A[0][0] = -vbn;
    A[0][1] = nx;
    A[0][2] = ny;
    A[0][3] = nz;
    A[0][4] = 0.;

    A[1][0] = ( -1. / pow(q[0],2) ) * ( q[1]*qn - gs*kq*nx );
    A[1][1] = ( 1. / q[0] ) * (q[1]*nx*(3.-GAMMA) + q[2]*ny + q[3]*nz ) - vbn;
    A[1][2] = ( 1. / q[0] ) * (q[1]*ny - gs*q[2]*nx );
    A[1][3] = ( 1. / q[0] ) * (q[1]*nz - gs*q[3]*nx );
    A[1][4] = gs * nx;

    A[2][0] = ( -1. / pow(q[0],2) ) * ( q[2]*qn - gs*kq*ny );
    A[2][1] = ( 1. / q[0] ) * (q[2]*nx - gs*q[1]*ny );
    A[2][2] = ( 1. / q[0] ) * (q[1] * nx + q[2]*ny*(3.-GAMMA) + q[3]*nz ) - vbn;
    A[2][3] = ( 1. / q[0] ) * (q[2]*nz - gs*q[3]*ny );
    A[2][4] = gs * ny;

    A[3][0] = ( -1. / pow(q[0],2) ) * ( q[3]*qn - gs*kq*nz );
    A[3][1] = ( 1. / q[0] ) * (q[3]*nx - gs*q[1]*nz );
    A[3][2] = ( 1. / q[0] ) * (q[3]*ny - gs*q[2]*nz );
    A[3][3] = ( 1. / q[0] ) * (q[1]*nx + q[2]*ny + q[3]*nz*(3.-GAMMA) ) - vbn;
    A[3][4] = gs * nz;

    A[4][0] = -( q[4]*qn / pow(q[0],2) ) * GAMMA + 2.*qn*gs*kq/pow(q[0],3);
    A[4][1] = ( q[4]/q[0] ) * GAMMA * nx - gs*q[1]/q[0];
    A[4][2] = ( q[4]/q[0] ) * GAMMA * ny - gs*q[2]/q[0];
    A[4][3] = ( q[4]/q[0] ) * GAMMA * nz - gs*q[3]/q[0];
    A[4][4] = ( qn / q[0] ) * GAMMA - vbn;
    
    return A;
}

void setStates(Vector<N_VAR>& primL, Vector<N_VAR>& consL, Vector<N_VAR>& primR, Vector<N_VAR>& consR, const Cell& LC, const Cell& RC, const CVector& disL, const CVector& disR)
{
    double k, ie;
    Vector2D<3,N_VAR> grad;
    
    minMod(LC.grad, RC.grad, grad);
    
    for (int i=0; i<N_VAR; ++i)
    {
        primL[i] = LC.prim[i] + dotP(grad[i], disL);
        primR[i] = RC.prim[i] + dotP(grad[i], disR);
    }
    
    //-----------------------------------------------------------------------
    
    consL[0] = primL[0];
    consL[1] = primL[0] * primL[1];
    consL[2] = primL[0] * primL[2];
    consL[3] = primL[0] * primL[3];

    k = 0.5 * ( pow(primL[1],2) + pow(primL[2],2) + pow(primL[3],2) );
    ie = primL[4] / ( (GAMMA-1.) * primL[0] );

    consL[4] = primL[0] * (k + ie);
    
    //-----------------------------------------------------------------------
    
    consR[0] = primR[0];
    consR[1] = primR[0] * primR[1];
    consR[2] = primR[0] * primR[2];
    consR[3] = primR[0] * primR[3];

    k = 0.5 * ( pow(primR[1],2) + pow(primR[2],2) + pow(primR[3],2) );
    ie = primR[4] / ( (GAMMA-1.) * primR[0] );

    consR[4] = primR[0] * (k + ie);
}

void roeflx(Face& face, Cell& LC, Cell& RC)
{
    double rho, u, v, w, H, a, k, qn, ql, qm;
    double RT;
    double drho, dp, dqn;
    double rhoL, uL, vL, wL, EL, kL, pL, qnL, qlL, qmL, HL, aL, ieL;
    double rhoR, uR, vR, wR, ER, kR, pR, qnR, qlR, qmR, HR, aR, ieR;
    double gamStar = GAMMA - 1.;
    double vb;
    double vbx;
    double vby;
    double vbz;
    double mg;
    double dql, dqm;
    double nx, ny, nz;
    double lx, ly, lz;
    double mx, my, mz;
    double vel;
    double asq;
    Vector<N_VAR> LdU;
    Vector<N_VAR> ws;
    Vector<N_VAR> diss;
    Vector<N_VAR> fluxL;
    Vector<N_VAR> fluxR;
    Vector<N_VAR> flux;
    Vector<N_VAR> primL;
    Vector<N_VAR> consL;
    Vector<N_VAR> primR;
    Vector<N_VAR> consR;
    CVector area;
    CVector disL;
    CVector disR;
    CVector n;
    Vector2D <N_VAR, N_VAR> R;
    Vector2D <N_VAR, N_VAR> L;
    Vector2D <N_VAR, N_VAR> Aroe;
    Vector2D <N_VAR, N_VAR> JL;
    Vector2D <N_VAR, N_VAR> JR;    
    BC bc;
    face_t boutype;    
    //--------------------------
    
    //Cell& LC = cell[ face.nei[0] ]; // left cell
    //Cell& RC = cell[ face.nei[1] ]; // right cell

    area = face.area;
    boutype = face.bouType;
    bc = RC.bc;
    
    n = norm(area);
    nx = n[0];
    ny = n[1];
    nz = n[2];
    
    mg = mag(area);

    vbx = face.vb[0];
    vby = face.vb[1];
    vbz = face.vb[2];

    // tangent vectors
    lx = ny;
    ly = -nx;
    lz = 0.;
    
    mx = -nz;
    my = 0.;
    mz = nx;
    
    disL = face.cnt - LC.cnt;
    disR = face.cnt - RC.cnt;
    
    setStates (primL, consL, primR, consR, LC, RC, disL, disR);

    // Left state
    rhoL = primL[0];
    uL   = primL[1];
    vL   = primL[2];
    wL   = primL[3];
    pL   = primL[4];
    ieL  = pL / ( gamStar * rhoL );
    kL   = 0.5 * (pow(uL,2) + pow(vL,2) + pow(wL,2));
    EL   = rhoL * (kL + ieL);
    aL   = sqrt(GAMMA*pL/rhoL);
    qnL  = uL*nx + vL*ny + wL*nz;
    qlL  = uL*lx + vL*ly + wL*lz;
    qmL  = uL*mx + vL*my + wL*mz;
    HL   = (EL + pL) / rhoL;

    // Right state
    rhoR = primR[0];
    uR   = primR[1];
    vR   = primR[2];
    wR   = primR[3];
    pR   = primR[4];
    ieR  = pR / ( gamStar * rhoR );
    kR   = 0.5 * (pow(uR,2) + pow(vR,2) + pow(wR,2));
    ER   = rhoR * (kR + ieR);
    
    aR   = sqrt(GAMMA*pR/rhoR);
    qnR  = uR*nx + vR*ny + wR*nz;
    qlR  = uR*lx + vR*ly + wR*lz;
    qmR  = uR*mx + vR*my + wR*mz;
    HR   = (ER + pR) / rhoR;

    vb = vbx*nx + vby*ny + vbz*nz;

    if (boutype == face_t::INTERIOR || bc == BC::DIRICHLET)
    {
        // Roe-averaged quantities
        RT  = sqrt(rhoR/rhoL);
        rho = RT * rhoL;
        u   = (uL + RT * uR) / (1. + RT);
        v   = (vL + RT * vR) / (1. + RT);
        w   = (wL + RT * wR) / (1. + RT);
        H   = (HL + RT * HR) / (1. + RT);
        k   = 0.5 * ( pow(u,2) + pow(v,2) + pow(w,2) );
        a   = sqrt( gamStar * (H - k) );
        qn  = u*nx + v*ny + w*nz;
        ql  = u*lx + v*ly + w*lz;
        qm  = u*mx + v*my + w*mz;
        asq = pow(a,2);

        // Wave strengths
        drho = rhoR - rhoL;
        dp   = pR - pL;
        dqn  = qnR - qnL;
        dql  = qlR - qlL;
        dqm  = qmR - qmL;

        LdU[0] = (dp - rho*a*dqn) / (2.*asq);
        LdU[1] = drho - dp/asq;
        LdU[2] = (dp + rho*a*dqn) / (2.*asq);
        LdU[3] = rho * dql;
        LdU[4] = rho * dqm;

        // Absolute values of the wave speeds
        ws[0] = fabs(qn - a  - vb);
        ws[1] = fabs(qn - vb);
        ws[2] = fabs(qn + a - vb);
        ws[3] = fabs(qn - vb);
        ws[4] = fabs(qn - vb);

        // Right eigenvectors
        R[0][0] = 1.;
        R[1][0] = u - a*nx;
        R[2][0] = v - a*ny;
        R[3][0] = w - a*nz;
        R[4][0] = H - a*qn;

        R[0][1] = 1.;
        R[1][1] = u;
        R[2][1] = v;
        R[3][1] = w;
        R[4][1] = k;

        R[0][2] = 1.;
        R[1][2] = u + a*nx;
        R[2][2] = v + a*ny;
        R[3][2] = w + a*nz;
        R[4][2] = H + a*qn;

        R[0][3] = 0.;
        R[1][3] = lx;
        R[2][3] = ly;
        R[3][3] = lz;
        R[4][3] = ql;

        R[0][4] = 0.;
        R[1][4] = mx;
        R[2][4] = my;
        R[3][4] = mz;
        R[4][4] = qm;

        // Dissipation
        for (int i=0; i<N_VAR; ++i)
        {
            diss[i] = 0.;
            for (int j=0; j<N_VAR; ++j)
            {
                diss[i] += ws[j] * LdU[j] * R[i][j];
            }
        }

        // Physical flux
        fluxL[0] = rhoL * (qnL - vb);
        fluxL[1] = rhoL * (qnL - vb) * uL + pL*nx;
        fluxL[2] = rhoL * (qnL - vb) * vL + pL*ny;
        fluxL[3] = rhoL * (qnL - vb) * wL + pL*nz;
        fluxL[4] = (qnL - vb) * EL + qnL * pL;

        fluxR[0] = rhoR * (qnR - vb);
        fluxR[1] = rhoR * (qnR - vb) * uR + pR*nx;
        fluxR[2] = rhoR * (qnR - vb) * vR + pR*ny;
        fluxR[3] = rhoR * (qnR - vb) * wR + pR*nz;
        fluxR[4] = (qnR - vb) * ER + qnR * pR;

        // Numerical flux
        for (int i=0; i<N_VAR; ++i)
        {
            flux[i] = 0.5 * (fluxL[i] + fluxR[i] - diss[i]) * mg;
        }

        vel = fabs(qn-vb) + a;
        
        //-------------------------------------------------------------
        
        // left eigenvectors (inv(R))
        L[0][0] = ( gamStar * k + a * qn) / (2. * asq);
        L[0][1] = (-gamStar * u - a * nx) / (2. * asq);
        L[0][2] = (-gamStar * v - a * ny) / (2. * asq);
        L[0][3] = (-gamStar * w - a * nz) / (2. * asq);
        L[0][4] = gamStar / (2. * asq);

        L[1][0] = (asq - gamStar * k) / asq;
        L[1][1] = (gamStar * u) / asq;
        L[1][2] = (gamStar * v) / asq;
        L[1][3] = (gamStar * w) / asq;
        L[1][4] = -gamStar / asq;

        L[2][0] = ( gamStar * k - a * qn) / (2. * asq);
        L[2][1] = (-gamStar * u + a * nx) / (2. * asq);
        L[2][2] = (-gamStar * v + a * ny) / (2. * asq);
        L[2][3] = (-gamStar * w + a * nz) / (2. * asq);
        L[2][4] = gamStar / (2. * asq);

        if ( fabs(nx) >= fabs(ny) && fabs(nx) >= fabs(nz) )
        {
            L[3][0] = (v - qn * ny) / nx;
            L[3][1] = ny;
            L[3][2] = (pow(ny,2) - 1.) / nx;
            L[3][3] = (ny * nz) / nx;
            L[3][4] = 0.;

            L[4][0] = (qn * nz - w) / nx;
            L[4][1] = -nz;
            L[4][2] = -(ny * nz) / nx;
            L[4][3] = (1. - pow(nz,2)) / nx;
            L[4][4] = 0.;
        }
        else if ( fabs(ny) >= fabs(nx) && fabs(ny) >= fabs(nz) )
        {
            L[3][0] = (qn * nx - u) / ny;
            L[3][1] = (1. - pow(nx,2)) / ny;
            L[3][2] = -nx;
            L[3][3] = -(nx * nz) / ny;
            L[3][4] = 0.;

            L[4][0] = (w - qn * nz) / ny;
            L[4][1] = (nx * nz) / ny;
            L[4][2] = nz;
            L[4][3] = (pow(nz,2) - 1.) / ny;
            L[4][4] = 0.;
        }
        else if ( fabs(nz) >= fabs(nx) && fabs(nz) >= fabs(ny) )
        {
            L[3][0] = (u - qn * nx) / nz;
            L[3][1] = (pow(nx,2) - 1.) / nz;
            L[3][2] = (nx * ny) / nz;
            L[3][3] = nx;
            L[3][4] = 0.;

            L[4][0] = (qn * ny - v) / nz;
            L[4][1] = -(nx * ny) / nz;
            L[4][2] = (1. - pow(ny,2)) / nz;
            L[4][3] = -ny;
            L[4][4] = 0.;
        }
        
        Aroe = mul5 ( mulI5 (R,ws), L); // Aroe = (R * ws) * L

        // left and right Jacobians
        JL = jacob(consL, n, vb);
        JR = jacob(consR, n, vb);
        /*JL = jacob(LC.cons, n, vb);
        JR = jacob(RC.cons, n, vb);*/

        // left and right matrices
        face.M[0] = add5(JL, Aroe); // AL = JL + Aroe
        face.M[1] = sub5(JR, Aroe); // AR = JR - Aroe
        
        face.M[0] = mulS5 (0.5*mg, face.M[0]);
        face.M[1] = mulS5 (0.5*mg, face.M[1]);
    }
    else if (bc == BC::SLIP_WALL)
    {
        flux[0] = rhoR * (qnR - vb);
        flux[1] = rhoR * (qnR - vb) * uR + pR*nx;
        flux[2] = rhoR * (qnR - vb) * vR + pR*ny;
        flux[3] = rhoR * (qnR - vb) * wR + pR*nz;
        flux[4] = (qnR - vb) * ER + qnR * pR;

        flux[0] *= mg;
        flux[1] *= mg;
        flux[2] *= mg;
        flux[3] *= mg;
        flux[4] *= mg;

        vel = fabs(qnR-vb) + aR;
        
        //-----------------------------------------------
        
        face.M[0] = jacob(consL, n, vb);
        face.M[1] = jacob(consR, n, vb);
        /*face.M[0] = jacob(LC.cons, n, vb);
        face.M[1] = jacob(RC.cons, n, vb);*/
        
        face.M[0] = mulS5 (mg, face.M[0]);
        face.M[1] = mulS5 (mg, face.M[1]);
        
        Vector2D <N_VAR,N_VAR> M;

        M[0][0] = 1.;
        M[0][1] = 0.;
        M[0][2] = 0.;
        M[0][3] = 0.;
        M[0][4] = 0.;

        M[1][0] = 0.;
        M[1][1] = 1. - 2. * pow(nx,2);
        M[1][2] = -2. * nx * ny;
        M[1][3] = -2. * nx * nz;
        M[1][4] = 0.;

        M[2][0] = 0.;
        M[2][1] = -2. * nx * ny;
        M[2][2] = 1. - 2. * pow(ny,2);
        M[2][3] = -2. * ny * nz;
        M[2][4] = 0.;

        M[3][0] = 0.;
        M[3][1] = -2. * nx * nz;
        M[3][2] = -2. * ny * nz;
        M[3][3] = 1. - 2. * pow(nz,2);
        M[3][4] = 0.;

        M[4][0] = 0.;
        M[4][1] = 0.;
        M[4][2] = 0.;
        M[4][3] = 0.;
        M[4][4] = 1.;
        
        face.M[0] = add5 ( mul5 (face.M[1], M), face.M[0] );
    }
    else if (bc == BC::EMPTY)
    {
        flux[0] = rhoR * (qnR - vb);
        flux[1] = rhoR * (qnR - vb) * uR + pR*nx;
        flux[2] = rhoR * (qnR - vb) * vR + pR*ny;
        flux[3] = rhoR * (qnR - vb) * wR + pR*nz;
        flux[4] = (qnR - vb) * ER + qnR * pR;

        flux[0] *= mg;
        flux[1] *= mg;
        flux[2] *= mg;
        flux[3] *= mg;
        flux[4] *= mg;

        vel = fabs(qnR-vb) + aR;
        
        //-----------------------------------------------
        
        face.M[0] = jacob(consL, n, vb);
        //face.M[0] = jacob(LC.cons, n, vb);
        face.M[0] = mulS5 (2.*mg, face.M[0]);
    }
    
    for (int i=0; i<N_VAR; ++i)
    {
        LC.R[i] -= flux[i];
        RC.R[i] += flux[i];
    }
    
    LC.sigma += mg * vel;
    RC.sigma += mg * vel;    
    
    LC.D = add5(LC.D, face.M[0]);
    RC.D = sub5(RC.D, face.M[1]);
}