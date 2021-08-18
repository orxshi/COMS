#include "AFT.h"

namespace AFT
{
    void createBoundaryElements (int c, const vector<Point>& points, const vector<Point>& points2, Grid& newGrid, const vector<Triangle>& triangles, int phys, int newGridId, vector <bool>& ptsExist, vector <int>& iGridPt)
    {
        int iP0, iP1, iP2, iP3;
        CVector na, iv;
        double dp;
        
        for (int t=0; t<triangles.size(); ++t)
        {
            Cell cll;
            
            cll.type = elmType_t::TRI;
            cll.phys = phys;
            cll.iBlank = iBlank_t::NA;
            cll.belonging = newGridId;
            cll.bc = BC::EMPTY;
            cll.vtx.reserve( static_cast<int>(nVertices_t::TRI) );
            cll.vtxBelo.resize (static_cast<int>(nVertices_t::TRI), newGridId);
            cll.ghost = true;
            //cll.newlyCreated = true;
            
            iP0 = triangles[t].p[0];
            iP1 = triangles[t].p[1];
            iP2 = triangles[t].p[2];
            
            ::Point pts[3];
            if (c < triangles.size())
            {
                pts[0].dim = points[iP0].dim;
                pts[1].dim = points[iP1].dim;
                pts[2].dim = points[iP2].dim;
            }
            else
            {   
                pts[0].dim = points2[iP0].dim;
                pts[1].dim = points2[iP1].dim;
                pts[2].dim = points2[iP2].dim;
            }
            
            // these calculations must be done before moving pts to vtx
            na = crossP( (pts[1].dim - pts[0].dim) , (pts[2].dim - pts[0].dim) );
            if (c < triangles.size()) { iv = points2[iP0].dim - pts[0].dim; }
            else                      { iv = points[iP0].dim - pts[0].dim; }
            dp = dotP( na , iv );            
            
            for (int i=0; i<triangles[t].p.size(); ++i)
            {
                int ip = triangles[t].p[i];
                
                if (c >= triangles.size()) { ip += points.size(); }
                
                if (ptsExist[ip])
                {
                    cll.vtx.push_back( iGridPt[ip] );
                }
                else
                {
                    ptsExist[ip] = true;                    
                    newGrid.pt.push_back( move(pts[i]) );
                    cll.vtx.push_back( newGrid.pt.size()-1 );
                    iGridPt[ip] = newGrid.pt.size()-1;
                }
            }

            if ((c < triangles.size() && dp < 0.) || (c >= triangles.size() && dp > 0.))
            {
                int tmp = cll.vtx[1];
                cll.vtx[1] = cll.vtx[2];
                cll.vtx[2] = tmp;
            }
            else if ( dp == 0. )
            {
                cout << "dp = 0." << endl;
                cout << "c = " << c << endl;
                cout << "triangles.size() = " << triangles.size() << endl;
                exit(-2);
            }
            
            // assign face
            Face nf;
            nf.vtx.reserve( 3 );
            nf.vtx.push_back( cll.vtx[0] );
            
            if (c < triangles.size())
            {
                nf.vtx.push_back( cll.vtx[2] );
                nf.vtx.push_back( cll.vtx[1] );
            }
            else
            {
                nf.vtx.push_back( cll.vtx[1] );
                nf.vtx.push_back( cll.vtx[2] );
            }
            
            nf.bouType = face_t::BOUNDARY;
            newGrid.face.push_back( nf );
            cll.face.push_back ( newGrid.face.size()-1 );
            newGrid.face.back().nei.push_back (c);
            
            newGrid.cell.push_back (move(cll));
            
            ++c;
        }
    }
    
    void createInternalFaceofCell (Grid& newGrid, int conf, Cell& cll, int t)
    {
        Face nf;
        nf.vtx.reserve( 4 );
        nf.bouType = face_t::INTERIOR;
        
        switch (conf)
        {
            case 0:
                nf.vtx.push_back( cll.vtx[0] );
                nf.vtx.push_back( cll.vtx[1] );
                nf.vtx.push_back( cll.vtx[4] );
                nf.vtx.push_back( cll.vtx[3] );
            break;
            
            case 1:
                nf.vtx.push_back( cll.vtx[2] );
                nf.vtx.push_back( cll.vtx[0] );
                nf.vtx.push_back( cll.vtx[3] );
                nf.vtx.push_back( cll.vtx[5] );
            break;
            
            case 2:
                nf.vtx.push_back( cll.vtx[1] );
                nf.vtx.push_back( cll.vtx[2] );
                nf.vtx.push_back( cll.vtx[5] );
                nf.vtx.push_back( cll.vtx[4] );
            break;
            
            default:
                cout << "undefined configuration in AFT::createInternalFaceofCell(...)" << endl;
                exit(-2);
            break;
        }
        
        int index;
        if (faceExists (nf, newGrid.face, newGrid.pt, index) == false)
        {
            newGrid.face.push_back( move(nf) );
            cll.face.push_back ( newGrid.face.size()-1 );
            newGrid.face.back().nei.push_back (t);
        }
        else
        {
            cll.face.push_back ( index );
            newGrid.face[index].nei.push_back (t);
        }
    }

    void createCells (double offsetZ, const vector<Point>& points, Grid& newGrid, const vector<Triangle>& triangles, int phys, int newGridId)
    {
        vector<Point> points2;
        points2.resize( points.size() );
        vector <bool> ptsExist (2*points.size(),false);// array.fill (false);
        vector <int> iGridPt (2*points.size(),-1);// array.fill (-1);

        // create offset points
        for (unsigned int i=0; i<points.size(); ++i)
        {
            points2[i] = points[i];
            points2[i].dim[2] = offsetZ;
        }
        
        newGrid.n_bou_elm = 2. * triangles.size();
        newGrid.totalNElms = newGrid.n_bou_elm + triangles.size();
        newGrid.n_in_elm = newGrid.totalNElms - newGrid.n_bou_elm;
        //newGrid.cell.resize( newGrid.totalNElms ); // better to reserve and push_back
        newGrid.cell.reserve (newGrid.totalNElms);

        cout << "creating boundary elements" << flush;
        
        createBoundaryElements (0, points, points2, newGrid, triangles, phys, newGridId, ptsExist, iGridPt);
        createBoundaryElements (triangles.size(), points, points2, newGrid, triangles, phys, newGridId, ptsExist, iGridPt);

        cout << "done!" << endl;

        // cells
        for (int t=newGrid.n_bou_elm; t<newGrid.totalNElms; ++t)
        {
            Cell cll;
        
            cll.type = elmType_t::PEN;
            cll.vtx.reserve( static_cast<int>(nVertices_t::PEN) );
            cll.vtxBelo.resize (static_cast<int>(nVertices_t::PEN), newGridId);
            cll.face.reserve( static_cast<int>(nFaces_t::PEN) );
            cll.iBlank = iBlank_t::FIELD;
            cll.belonging = newGridId;
            cll.ghost = false;
            //cll.newlyCreated = true;
            
            const Cell& noCell = newGrid.cell[t-newGrid.n_bou_elm];
            const Cell& oCell  = newGrid.cell[t-triangles.size()];

            // assign vertices
            cll.vtx.push_back( noCell.vtx[0] );
            cll.vtx.push_back( noCell.vtx[1] );
            cll.vtx.push_back( noCell.vtx[2] );
            cll.vtx.push_back( oCell.vtx[0] );
            cll.vtx.push_back( oCell.vtx[1] );
            cll.vtx.push_back( oCell.vtx[2] );
            
            // assign face from non-offset boundary element
            cll.face.push_back ( noCell.face[0] );
            newGrid.face[ noCell.face[0] ].nei.push_back (t);
            
            // assign face from offset boundary element
            cll.face.push_back ( oCell.face[0] );
            newGrid.face[ oCell.face[0] ].nei.push_back (t);
            
            createInternalFaceofCell (newGrid, 0, cll, t);
            createInternalFaceofCell (newGrid, 1, cll, t);
            createInternalFaceofCell (newGrid, 2, cll, t);
            
            newGrid.cell.push_back (move(cll));
        }
        
        for (Cell& cll: newGrid.cell)
        {
            cll.set_centroid(newGrid.pt);
        }
        
        for (Face& f: newGrid.face)
        {
            f.set_area (newGrid.pt);
            f.set_centroid (newGrid.pt);
            
            if (f.nei.size() == 2)
            {
                newGrid.cell[ f.nei[0] ].nei.push_back( f.nei[1] );
                newGrid.cell[ f.nei[1] ].nei.push_back( f.nei[0] );
            }
            
            if (f.bouType == face_t::BOUNDARY)
            {
                if (f.nei.size() != 2)
                {
                    cout << "f.nei.size() != 2 in createCells(...)" << endl;
                    exit(-2);
                }
                
                int tmp = f.nei[0];
                f.nei[0] = f.nei[1];
                f.nei[1] = tmp;
            }
        }
        
        newGrid.set_elmVolumes();
    }
    
    bool faceExists (const Face& nf, const vector<Face>& face, const vector<::Point>& point, int& index)
    {
        int counter;
        
        for (unsigned int iF=0; iF<face.size(); ++iF)
        {
            const Face& f = face[iF];
            
            counter = 0;
            for (int iV=0; iV<f.vtx.size(); ++iV)
            {
                for (int iNF=0; iNF<nf.vtx.size(); ++iNF)
                {
                    if (f.vtx[iV] == nf.vtx[iNF])
                    {
                        ++counter;
                    }                    
                }
            }
            
            if (counter == f.vtx.size())
            {
                index = iF;
                return true;
            }
        }
        
        return false;
    }
}
