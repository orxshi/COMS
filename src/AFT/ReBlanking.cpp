#include "AFT.h"

void Grid::fieldToFringe (int crt)
{
    int counter;

    for (Cell& c: cell)
    {
        if (c.iBlank == iBlank_t::FIELD)
        {
            counter = 0;
            for (const int nei: c.nei)
            {
                Cell& n = cell[nei];
                
                if (n.iBlank == iBlank_t::FRINGE)
                {
                    ++counter;
                    if (counter == crt)
                    {
                        for (int i=0; i<c.receiver.size(); ++i)
                        {
                            c.receiver[i]->iBlank = iBlank_t::FIELD;
                            c.receiver[i]->donor = NULL;
                            c.receiver[i]->receiver.push_back (&c);
                        }
                        

                        c.iBlank = iBlank_t::FRINGE;
                        c.donor = c.receiver[0]; // may do better choice
                        
                        for (int i=0; i<c.receiver.size(); ++i)
                        {
                            c.receiver[i] = NULL;
                        }
                        c.receiver.clear();
                        
                        break;
                    }
                }
            }
        }
    }
}
    
void Grid::fringeToField (int crt)
{
    int counter;    

    for (Cell& c: cell)
    {
        if (c.iBlank == iBlank_t::FRINGE)
        {
            counter = 0;            
            for (const int nei: c.nei)
            {
                Cell& n = cell[nei];
                
                if (n.iBlank == iBlank_t::FIELD)
                {
                    ++counter;
                    if (counter == crt)
                    {
                        //c.donor->iBlank = iBlank_t::FRINGE;
                        //c.donor->donor = &c;
                        
                        for (int i=0; i<c.donor->receiver.size(); ++i)
                        {
                            if (c.donor->receiver[i] == &c)
                            {
                                c.donor->receiver[i] = NULL;
                                c.donor->receiver.erase (c.donor->receiver.begin()+i);
                                break;
                            }
                        }

                        c.iBlank = iBlank_t::FIELD;
                        c.donor = NULL;

                        break;
                    }
                }
            }
        }
    }
}

/*void Grid::blankWithPhys (int phys)
{
    for (int c=0; c<n_bou_elm; ++c)
    {
        Cell& cll = cell[c];
        
        if (cll.phys == phys)
        {
            cll.iBlank = iBlank_t::FRINGE;
            cell[cll.nei[0]].iBlank = iBlank_t::FRINGE;
            // but who is donor???
        }
    }
}*/