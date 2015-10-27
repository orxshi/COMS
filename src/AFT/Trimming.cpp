#include "AFT.h"

void Grid::trimWhoHasFringeNeighbor()
{
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        Cell& cll = cell[c];
        
        if (cll.iBlank == iBlank_t::FIELD)
        {
            for (const int n: cll.nei)
            {
                const Cell& nei = cell[ n ];
                
                if (nei.iBlank == iBlank_t::FRINGE)
                {
                    cll.trim = true;
                    ++cll.nTrims;
                    break;
                }
            }
        }
    }
}

void Grid::trimWhoHasTrimNeighbor (int threshold)
{
    int counter;
    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        Cell& cll = cell[c];
        
        if (cll.iBlank == iBlank_t::FIELD)
        {
            counter = 0;
            
            for (const int n: cll.nei)
            {
                const Cell& nei = cell[n];
                
                if (nei.trim == true && nei.iBlank == iBlank_t::FIELD && nei.nTrims == 1)
                {
                    ++counter;
                    
                    if (counter == threshold)
                    {
                        cll.trim = true;
                        ++cll.nTrims;
                        break;
                    }
                }
            }
        }
    }
}

void Grid::trimToUntrim (int crt)
{
    int counter;

    for (int ic=n_bou_elm; ic<cell.size(); ++ic)
    {
        Cell& c = cell[ic];
        
        if (c.trim == true)
        {
            counter = 0;
            for (const int nei: c.nei)
            {
                Cell& n = cell[nei];
                
                if (n.iBlank == iBlank_t::FIELD && n.trim == false)
                {
                    ++counter;
                    if (counter == crt)
                    {
                        c.trim = false;
                        break;
                    }
                }
            }
        }
    }
}

void Grid::nontrimToTrim (int crt)
{
    int counter;

    for (int ic=n_bou_elm; ic<cell.size(); ++ic)
    {
        Cell& c = cell[ic];
        
        if (c.trim == false)
        {
            counter = 0;
            for (const int nei: c.nei)
            {
                Cell& n = cell[nei];
                
                if (n.trim == true)
                {
                    ++counter;
                    if (counter == crt)
                    {
                        c.trim = true;
                        break;
                    }
                }
            }
        }
    }
}

/*void aft::TrimWhoHasTrimNeighborWithFringeNeighbour (Grid& grid)
{
    int nei, neii;

    for (unsigned int e=grid.n_bou_elm; e<grid.n_elm; ++e)
    {
        if (grid.elm[e].iBlank == 1)
        {
            for (unsigned int n=0; n<grid.elm[e].n_faces; ++n)
            {
                nei = grid.elm[e].neigh[n];

                if (grid.elm[nei].trim == 1)
                {
                    for (unsigned int nn=0; nn<grid.elm[nei].n_faces; ++nn)
                    {
                        neii = grid.elm[nei].neigh[nn];

                        if (grid.elm[neii].iBlank == 0)
                        {
                            grid.elm[e].trim = 1;
                            break;
                        }
                    }
                }

                if (grid.elm[e].trim == 1)
                {
                    break;
                }
            }
        }
    }
}*/

