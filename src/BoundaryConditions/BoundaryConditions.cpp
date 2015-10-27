#include "../Grid/Grid.h"

void empty (Cell& cll, const vector<Cell>& cell)
{
    cll.prim = cell[ cll.nei[0] ].prim;
}

void slipWall (Cell& cll, const vector<Face>& face, const vector<Cell>& cell)
{
    CVector vel, normVel, tangVel;
    
    const Cell& nei = cell[ cll.nei[0] ];
    const Face& f = face[ cll.face[0] ];
    CVector n = norm(f.area);
    
    cll.prim[0] = nei.prim[0];
    cll.prim[4] = nei.prim[4];

    vel[0] = nei.prim[1];
    vel[1] = nei.prim[2];
    vel[2] = nei.prim[3];
    
    normVel = dotP(vel,n) * n; // normal velocity
    tangVel = vel - normVel; // tangential velocity

    cll.prim[1] = tangVel[0] + f.vb[0];
    cll.prim[2] = tangVel[1] + f.vb[1];
    cll.prim[3] = tangVel[2] + f.vb[2];
}

void Grid::apply_BCs()
{
    for (int c=0; c<n_bou_elm; ++c)
    {
        switch (cell[c].bc)
        {
            case BC::EMPTY:
                empty( cell[c],cell );
                break;
            case BC::SLIP_WALL:
                slipWall( cell[c],face,cell );
                break;
        }
        
        cell[c].prim_to_cons();
    }
}
void Grid::set_BCs()
{
    for (Cell& cll: cell)
    {
        cll.bc = BC::NA;
    }
    
    for (int c=0; c<n_bou_elm; ++c)
    {
        for (int i=0; i<phys_count; ++i)
        {
            if (cell[c].phys == phys[i])
            {
                cell[c].bc = static_cast<BC> (bc[i]);
                break;
            }
        }
    }
    
    
}
