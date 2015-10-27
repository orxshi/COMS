#include "Grid.h"

void Grid::set_faceVertices (Face& face, const Cell& elm, const int index)
{
    switch (elm.type)
    {
        case elmType_t::HEX:

            switch (index)
            {
                
                
                case 0:
                    face.vtx.push_back ( elm.vtx[0] );
                    face.vtx.push_back ( elm.vtx[3] );
                    face.vtx.push_back ( elm.vtx[2] );
                    face.vtx.push_back ( elm.vtx[1] );
                    break;

                case 1:
                    face.vtx.push_back ( elm.vtx[0] );
                    face.vtx.push_back ( elm.vtx[1] );
                    face.vtx.push_back ( elm.vtx[5] );
                    face.vtx.push_back ( elm.vtx[4] );
                    break;

                case 2:
                    face.vtx.push_back ( elm.vtx[1] );
                    face.vtx.push_back ( elm.vtx[2] );
                    face.vtx.push_back ( elm.vtx[6] );
                    face.vtx.push_back ( elm.vtx[5] );
                    break;

                case 3:
                    face.vtx.push_back ( elm.vtx[2] );
                    face.vtx.push_back ( elm.vtx[3] );
                    face.vtx.push_back ( elm.vtx[7] );
                    face.vtx.push_back ( elm.vtx[6] );
                    break;

                case 4:
                    face.vtx.push_back ( elm.vtx[0] ); 
                    face.vtx.push_back ( elm.vtx[4] ); 
                    face.vtx.push_back ( elm.vtx[7] ); 
                    face.vtx.push_back ( elm.vtx[3] ); 
                    break;

                case 5:
                    face.vtx.push_back ( elm.vtx[4] ); 
                    face.vtx.push_back ( elm.vtx[5] ); 
                    face.vtx.push_back ( elm.vtx[6] ); 
                    face.vtx.push_back ( elm.vtx[7] ); //face.vtxIndex.push_back (elm.vtxIndex[7]);
                    break;
            }
            break;

        case elmType_t::TET:

            switch (index)
            {
                
                
                case 0:
                    face.vtx.push_back ( elm.vtx[0] ); 
                    face.vtx.push_back ( elm.vtx[2] ); 
                    face.vtx.push_back ( elm.vtx[1] ); 
                    break;

                case 1:
                    face.vtx.push_back ( elm.vtx[0] ); 
                    face.vtx.push_back ( elm.vtx[1] ); 
                    face.vtx.push_back ( elm.vtx[3] ); 
                break;

                case 2:
                    face.vtx.push_back ( elm.vtx[1] ); 
                    face.vtx.push_back ( elm.vtx[2] ); 
                    face.vtx.push_back ( elm.vtx[3] ); 
                break;

                case 3:
                    face.vtx.push_back ( elm.vtx[2] ); 
                    face.vtx.push_back ( elm.vtx[0] ); 
                    face.vtx.push_back ( elm.vtx[3] ); 
                break;
            }
            break;

        case elmType_t::PEN:

        switch(index)
        {
            
            
            case 0:
                face.vtx.push_back ( elm.vtx[0] ); 
                face.vtx.push_back ( elm.vtx[1] ); 
                face.vtx.push_back ( elm.vtx[4] ); 
                face.vtx.push_back ( elm.vtx[3] ); 
                break;

            case 1:
                face.vtx.push_back ( elm.vtx[1] ); 
                face.vtx.push_back ( elm.vtx[2] ); 
                face.vtx.push_back ( elm.vtx[5] ); 
                face.vtx.push_back ( elm.vtx[4] ); 
                break;

            case 2:
                face.vtx.push_back ( elm.vtx[2] ); 
                face.vtx.push_back ( elm.vtx[0] ); 
                face.vtx.push_back ( elm.vtx[3] ); 
                face.vtx.push_back ( elm.vtx[5] ); 
                break;

            case 3:
                face.vtx.push_back ( elm.vtx[0] ); 
                face.vtx.push_back ( elm.vtx[2] );
                face.vtx.push_back ( elm.vtx[1] );
                break;

            case 4:
                face.vtx.push_back ( elm.vtx[3] ); 
                face.vtx.push_back ( elm.vtx[4] ); 
                face.vtx.push_back ( elm.vtx[5] ); 
                break;
        }
        break;
        
        default:
            cout << "Cell type = " << static_cast<int>(elm.type) << endl;
            exit(-2);
            break;
    }
}

void Grid::set_connectivity()
{
    /**
    face.LR
    face.area
    face.centroid
    face.boutype

    elm.neigh
    elm.face
    elm.vol

    are set here.
    */

    vector <CVector> vtx;
    vector <faceNumCheck> fnc (totalNElms);
    unsigned int sig;
    double tmpd;

    for (unsigned int e=n_bou_elm; e<cell.size(); ++e)
    {
        cell[e].face.reserve (cell[e].nFaces);
        
        for (unsigned int j=0; j<cell[e].nFaces; ++j)
        {
            sig = 0;

            Face tmpFace;
            
            if (cell[e].type == elmType_t::HEX)
            {
                tmpFace.vtx.reserve (6*4);                            
            }
            else if (cell[e].type == elmType_t::TET)
            {
                tmpFace.vtx.reserve (4*3);
            }
            else if (cell[e].type == elmType_t::PEN)
            {
                tmpFace.vtx.reserve (5*4);
            }
            
            set_faceVertices (tmpFace, cell[e], j);

            for (unsigned int nei: bt[tmpFace.vtx[0]].a)
            {
                //if (nei >= n_bou_elm)
                //{
                    if (nei == e)
                    {
                        continue;
                    }
                //}
                
                for (auto n=tmpFace.vtx.begin()+1; n<tmpFace.vtx.end(); ++n)
                {
                    sig += bt[*n].search(nei);
                }

                if (sig == tmpFace.vtx.size()-1 && fnc[e].nmap[nei] != 1)
                {
                    fnc[e].nmap[nei] = 1;
                    fnc[nei].nmap[e] = 1;
                    
                    cell[e]  .nei.push_back (nei);
                    cell[nei].nei.push_back (e);

                    tmpFace.nei.push_back (e);
                    tmpFace.nei.push_back (nei);

                    tmpFace.set_area (pt);
                    tmpFace.set_centroid(pt);
                    
                    if (nei < n_bou_elm)
                    {
                        tmpFace.bouType = face_t::BOUNDARY;
                    }
                    else
                    {
                        tmpFace.bouType = face_t::INTERIOR;
                    }

                    // set volumes
                    tmpd = tmpFace.cnt[0] * tmpFace.area[0];
                    cell[e]  .vol += tmpd;
                    cell[nei].vol -= tmpd;
                    //tmpFace.nei[0].get().vol += tmpd;
                    //tmpFace.nei[1].get().vol -= tmpd;
                    
                    face.push_back( tmpFace );
                    
                    cell[e]  .face.push_back( face.size()-1 );
                    cell[nei].face.push_back( face.size()-1 );

                    break;
                }

                sig = 0;
            }
        }
    }
    
    for (int c=0; c<n_bou_elm; ++c)
    {
        if (cell[c].nei.size() != 1)
        {
            cout << "!!! Error: cell[c].nei.size() != 1 in set_connectivity()" << endl;
            exit(-2);
        }
    }
}

void Grid::set_elmVolumes ()
{
    double tmpd;
    
    for (const Face& f: face)
    {
        tmpd = f.cnt[1] * f.area[1];
        
        if (f.nei.size() == 1)
        {
            cell[f.nei[0]].vol += tmpd;
        }
        else if (f.nei.size() == 2)
        {
            cell[f.nei[0]].vol += tmpd;
            cell[f.nei[1]].vol -= tmpd;
        }
        else
        {
            cout << "!!! Error: f.nei.size() != 1 or 2" << endl;
            exit(-2);
        }
    }
}

void Grid::set_elmCentroids ()
{   
    for (Cell& c: cell)
    {
        c.set_centroid(pt);
    }
}

void Grid::set_grid ()
{
    set_connectivity();
    bt.clear();
    bt.shrink_to_fit();
    set_elmCentroids();
    leastSquaresCoeffs();
}

