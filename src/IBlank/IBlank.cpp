#include "../Grid/Grid.h"
#include "IBlank.h"

void Grid::CellADT::build (const Grid& gr)
{
    points.resize (gr.n_in_elm);
    
    for (unsigned int c=0; c<points.size(); ++c)
    {
        for (int i=0; i<ADT_DIM; ++i)
        {
            points[c].dim[i*2]   = BIG_POS_NUM;
            points[c].dim[i*2+1] = BIG_NEG_NUM;
        }

        for (const int ip: gr.cell[c+gr.n_bou_elm].vtx)
        {
            const Point& p = gr.pt[ip];
            
            for (int i=0; i<ADT_DIM; ++i)
            {
                points[c].dim[i*2]   = min (p.dim[i], points[c].dim[i*2]);
                points[c].dim[i*2+1] = max (p.dim[i], points[c].dim[i*2+1]);
            }

            points[c].vertices.push_back (p.dim);
        }

        points[c].idx = c+gr.n_bou_elm;
    }
    
    ADT::build();
}

void Grid::setWallDistance (int phys)
{
    for (int c=0; c<cell.size(); ++c)
    {
        double d = BIG_POS_NUM;
        
        for (int g=0; g<n_bou_elm; ++g)
        {
            if (cell[g].phys == phys)
            {
                CVector tmpv = cell[g].cnt - cell[c].cnt;

                d = min ( d, mag (tmpv) );
            }
        }

        cell[c].wallDistance = d;
    }
}

Iblank::Iblank ()
{
    // default
    cellCriter = cellCriter_t::WALL;

    ifstream in;
    in.open("iblank.dat");
    
    if (in.is_open())
    {
        string tmps;
        int tmpi;
        in >> tmps; in >> tmpi;
        
        cellCriter = static_cast<cellCriter_t> (tmpi);
    }
    
    in.close();
}

void Iblank::identify (Grid& grAct, Grid& grPas)
{
    function<int(Cell&)> getIndex = [&] (Cell& cll)
    {
        ADT::ADTPoint vec;
        CVector cnt;
        
        cnt = cll.cnt;
        
        for (int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = cnt[i];
            vec.dim[i*2+1] = cnt[i];
        }

        /*vec.dim[0] = cnt[0];
        vec.dim[2] = cnt[1];
        vec.dim[4] = cnt[2];

        vec.dim[1] = cnt[0];
        vec.dim[3] = cnt[1];
        vec.dim[5] = cnt[2];*/        
        
        return grPas.cellADT.search (vec);
    };

    for (int c=grAct.n_bou_elm; c<grAct.cell.size(); ++c)
    {
        int index;        
        Cell& cll = grAct.cell[c];
        
        if (cll.iBlank == iBlank_t::UNDEFINED)
        {            
            index = getIndex (cll);
            
            if (index == -1)
            {
                bool insideHole;
                
                if (grPas.nHoles > 0)
                {
                    for (int i=0; i<grPas.nHoles; ++i)
                    {
                        for (int j=0; j<N_DIM; ++j)
                        {
                            if (cll.cnt[j] >= grPas.holes[i].min[j] && cll.cnt[j] <= grPas.holes[i].max[j])
                            {
                                insideHole = true;
                            }
                            else
                            {
                                insideHole = false;
                                break;
                            }
                        }
                    }                    
                }
                
                if (grPas.nHoles > 0 && insideHole)
                {
                    cll.iBlank = iBlank_t::HOLE;                    
                }
                else
                {
                    cll.iBlank = iBlank_t::FIELD;
                }
            }
            else
            {
                if (grPas.cell[index].iBlank == iBlank_t::FIELD)
                {
                    cll.iBlank = iBlank_t::FRINGE;
                    cll.donor = &grPas.cell[index];
                    grPas.cell[index].receiver.push_back (&cll);
                }
                else if (grPas.cell[index].iBlank == iBlank_t::FRINGE)
                {
                    cll.iBlank = iBlank_t::FIELD;
                }
                else if (grPas.cell[index].iBlank == iBlank_t::UNDEFINED)
                {
                    double var1;
                    double var2;
                    
                    if (cellCriter == cellCriter_t::WALL)
                    {
                        var1 = grPas.cell[index].wallDistance;
                        var2 = cll.wallDistance;
                    }
                    else if (cellCriter == cellCriter_t::SIZE)
                    {
                        var1 = grPas.cell[index].vol;
                        var2 = cll.vol;
                    }
                    
                    if (var1 < var2)
                    {
                        cll.iBlank = iBlank_t::FRINGE;
                        grPas.cell[index].iBlank = iBlank_t::FIELD;
                        cll.donor = &grPas.cell[index];
                        grPas.cell[index].receiver.push_back (&cll);
                    }
                    else
                    {
                        cll.iBlank = iBlank_t::FIELD;
                        grPas.cell[index].iBlank = iBlank_t::FRINGE;
                        grPas.cell[index].donor = &cll;
                        cll.receiver.push_back (&grPas.cell[index]);
                    }
                }
                else
                {
                    cout << "undefined behavior in Grid::identifyIBlank(...)" << endl;
                    exit(-2);
                }
            }
        }
    }
    
    for (int c=0; c<grAct.n_bou_elm; ++c)
    {
        int index;
        Cell& cll = grAct.cell[c];
        
        if (cll.iBlank == iBlank_t::UNDEFINED)
        {
            if (cll.fringeBou == fringeBou_t::YES)
            {
                index = getIndex (cll);

                if (index == -1)
                {
                    cll.iBlank == iBlank_t::NA;
                }
                else
                {
                    if (grPas.cell[index].iBlank == iBlank_t::FIELD)
                    {
                        cll.iBlank = iBlank_t::FRINGE;
                        cll.donor = &grPas.cell[index];
                        grPas.cell[index].receiver.push_back (&cll);
                    }
                    else if (grPas.cell[index].iBlank == iBlank_t::FRINGE)
                    {                        
                        cll.iBlank = iBlank_t::FRINGE;
                        grPas.cell[index].iBlank = iBlank_t::FIELD;
                        
                        for (int r=0; r<grPas.cell[index].donor->receiver.size(); ++r)
                        {
                            if (grPas.cell[index].donor->receiver[r] == &grPas.cell[index])
                            {
                                grPas.cell[index].donor->receiver[r] = NULL;
                                grPas.cell[index].donor->receiver.erase (grPas.cell[index].donor->receiver.begin() + r);
                                break;
                            }
                        }
                        
                        grPas.cell[index].donor = NULL;
                        grPas.cell[index].receiver.push_back (&cll);
                        cll.donor = &grPas.cell[index];
                        //cout << "unset situation in Grid::identifyIBlank(...)" << endl;
                        //exit(-2);
                    }
                    else if (grPas.cell[index].iBlank == iBlank_t::UNDEFINED)
                    {
                        cll.iBlank = iBlank_t::FRINGE;
                        cll.donor = &grPas.cell[index];
                        grPas.cell[index].receiver.push_back (&cll);
                        grPas.cell[index].iBlank = iBlank_t::FIELD;
                    }
                    else
                    {
                        cout << "second undefined behavior in Grid::identifyIBlank(...)" << endl;
                        cout << "iBlank = " << static_cast<int>(grPas.cell[index].iBlank) << endl;
                        cout << "bc = " << static_cast<int>(cll.bc) << endl;
                        cout << "cell[c].cnt[0] = " << cll.cnt[0] << endl;
                        cout << "cell[c].cnt[1] = " << cll.cnt[1] << endl;
                        cout << "cell[c].cnt[2] = " << cll.cnt[2] << endl;
                        cout << "cell[c].belonging = " << cll.belonging << endl;
                        cout << "c = " << c << endl;
                        exit(-2);
                    }
                }
            }
            else if (cll.fringeBou == fringeBou_t::NO)
            {
                cll.iBlank == iBlank_t::NA;
            }
            else if (cll.fringeBou == fringeBou_t::UNDEFINED)
            {
                cout << "undefined fringeBou_t in Grid::identifyIBlank(...)" << endl;
                exit(-2);
            }
        }
    }
}

bool Grid::CellADT::compareFunction (const Node* node, const ADTPoint& targetPoint)
{
    bool inside;
    CVector tempVector;
    unsigned int size = node->p->vertices.size();
    double frac[size];
    vector<CVector> xc(size);

    for (unsigned int v=0; v<size; ++v)
    {
        xc[v] = node->p->vertices[v];
    }
    
    for (int i=0; i<ADT_DIM; ++i)
    {
        tempVector[i] = targetPoint.dim[i*2];
    }
    
    /*tempVector[0] = targetPoint.dim[0];
    tempVector[1] = targetPoint.dim[2];
    tempVector[2] = targetPoint.dim[4];*/
    
    osInterpolants (xc, tempVector, size, frac);

    inside = true;
    for (unsigned int v=0; v<size; ++v)
    {
        if (frac[v]<=0 || frac[v]>=1)
        {
            inside = false;
            break;
        }
    }

    return inside;
}

bool Grid::CellADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
{
    bool insideCube = true;

    for (int d=0; d<ADT_DIM; ++d)
    {
        if (!(node->p->dim[d*2] <= targetPoint.dim[d*2]) || !(node->p->dim[d*2+1] >= targetPoint.dim[d*2+1]))
        {
            insideCube = false;
            break;
        }
    }

    return insideCube;
}

void Grid::interpolate()
{
    for (Cell& cll: cell) { cll.interpolate(); }
    apply_BCs();
}