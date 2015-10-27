#include "Cell.h"

Cell::Cell()
{
    vol    = 0.;
    iBlank = iBlank_t::UNDEFINED;
    trim   = false;
    newlyCreated = false;
    type = elmType_t::UNDEFINED;
    bc = BC::UNDEFINED;
    donor = NULL;
    //receiver = NULL;
    for (Vector<N_DIM>& i: grad)
    {
        i.fill (0.);
    }
    fringeBou = fringeBou_t::UNDEFINED;
    nTrims = 0;
    nVertices = 0;
    nFaces = 0;
    
    dQ[0] = 0.;
    dQ[1] = 0.;
    dQ[2] = 0.;
    dQ[3] = 0.;
    dQ[4] = 0.;
}

Cell::Cell (const Cell& other)
{
    phys = other.phys;
    donor = NULL;
    
    for (int i=0; i<receiver.size(); ++i)
    {
        receiver[i] = NULL;
    }
    receiver.clear();
    
    iBlank = other.iBlank;
    belonging = other.belonging;    
    nei = other.nei;
    Mach = other.Mach;
    sigma = other.sigma;
    wallDistance = other.wallDistance;
    vol = other.vol;
    r_11 = other.r_11;
    r_12 = other.r_12;
    r_22 = other.r_22;
    r_23 = other.r_23;
    r_33 = other.r_33;
    R = other.R;
    dQ = other.dQ;
    old_dQ = other.old_dQ;
    prim = other.prim;
    cons = other.cons;
    old_cons = other.old_cons;
    oldold_cons = other.oldold_cons;
    resInner = other.resInner;
    resOuter = other.resOuter;
    D = other.D;
    trim = other.trim;
    newlyCreated = other.newlyCreated;
    ghost = other.ghost;
    cnt = other.cnt;
    grad = other.grad;
    emin = other.emin;
    emax = other.emax;
    face = other.face;
    type = other.type;
    vtx = other.vtx;
    vtxBelo = other.vtxBelo;
    bc = other.bc;
    fringeBou = other.fringeBou;
    nTrims = other.nTrims;
    nVertices = other.nVertices;
    nFaces = other.nFaces;
}

Cell& Cell::operator= (const Cell& other)
{
    phys = other.phys;
    donor = NULL;
    
    for (int i=0; i<receiver.size(); ++i)
    {
        receiver[i] = NULL;
    }
    receiver.clear();

    iBlank = other.iBlank;
    belonging = other.belonging;    
    nei = other.nei;
    Mach = other.Mach;
    sigma = other.sigma;
    wallDistance = other.wallDistance;
    vol = other.vol;
    r_11 = other.r_11;
    r_12 = other.r_12;
    r_22 = other.r_22;
    r_23 = other.r_23;
    r_33 = other.r_33;
    R = other.R;
    dQ = other.dQ;
    old_dQ = other.old_dQ;
    prim = other.prim;
    cons = other.cons;
    old_cons = other.old_cons;
    oldold_cons = other.oldold_cons;
    resInner = other.resInner;
    resOuter = other.resOuter;
    D = other.D;
    trim = other.trim;
    newlyCreated = other.newlyCreated;
    ghost = other.ghost;
    cnt = other.cnt;
    grad = other.grad;
    emin = other.emin;
    emax = other.emax;
    face = other.face;    
    type = other.type;
    vtx = other.vtx;
    vtxBelo = other.vtxBelo;
    bc = other.bc;
    fringeBou = other.fringeBou;
    nTrims = other.nTrims;
    nVertices = other.nVertices;
    nFaces = other.nFaces;
}

Cell::Cell (Cell&& other) :
phys (move(other.phys)),
iBlank (move(other.iBlank)),
belonging(move(other.belonging)),
nei(move(other.nei)),
Mach(move(other.Mach)),
sigma(move(other.sigma)),
wallDistance(move(other.wallDistance)),
vol(move(other.vol)),
r_11(move(other.r_11)),
r_12(move(other.r_12)),
r_22(move(other.r_22)),
r_23(move(other.r_23)),
r_33(move(other.r_33)),
R(move(other.R)),
dQ(move(other.dQ)),
old_dQ(move(other.old_dQ)),
prim(move(other.prim)),
cons(move(other.cons)),
old_cons(move(other.old_cons)),
oldold_cons(move(other.oldold_cons)),
resInner(move(other.resInner)),
resOuter(move(other.resOuter)),
D(move(other.D)),
trim(move(other.trim)),
newlyCreated(move(other.newlyCreated)),
ghost(move(other.ghost)),
cnt(move(other.cnt)),
grad(move(other.grad)),
emin(move(other.emin)),
emax(move(other.emax)),   
face(move(other.face)),
type(move(other.type)),
vtx(move(other.vtx)),
vtxBelo(move(other.vtxBelo)),
bc(move(other.bc)),
fringeBou(move(other.fringeBou)),
nTrims(move(other.nTrims)),
nVertices(move(other.nVertices)),
nFaces(move(other.nFaces))
{    
    donor = NULL;
    
    for (int i=0; i<receiver.size(); ++i)
    {
        receiver[i] = NULL;
    }
    receiver.clear();
}

void Cell::set_centroid (const vector<Point>& pt)
{
    for (unsigned int i=0; i<cnt.size(); ++i)
    {
        cnt[i] = 0.;
    }
    
    for (const int p: vtx)
    {
        cnt += pt[p].dim;
    }

    cnt /= vtx.size();
}

void Cell::cons_to_prim()
{
    double k,ie;    
    
    prim[0] = cons[0];
    prim[1] = cons[1] / prim[0];
    prim[2] = cons[2] / prim[0];
    prim[3] = cons[3] / prim[0];

    k = 0.5 * ( pow(prim[1],2) + pow(prim[2],2) + pow(prim[3],2) );
    ie = cons[4] / prim[0] - k;

    prim[4] = prim[0] * (GAMMA - 1.) * ie;
}

void Cell::prim_to_cons()
{
    double k,ie;
    
    cons[0] = prim[0];
    cons[1] = prim[0] * prim[1];
    cons[2] = prim[0] * prim[2];
    cons[3] = prim[0] * prim[3];

    k = 0.5 * ( pow(prim[1],2) + pow(prim[2],2) + pow(prim[3],2) );
    ie = prim[4] / ( (GAMMA-1.) * prim[0] );

    cons[4] = prim[0] * (k + ie);
}

void Cell::interpolate()
{
    CVector dis;        
    
    if ( iBlank == iBlank_t::FRINGE )
    {
        for (int i=0; i<N_VAR; ++i)
        {
            dis = cnt - (*donor).cnt;

            prim[i] = donor->prim[i] + dotP(donor->grad[i], dis);
        }

        prim_to_cons();
    }
}