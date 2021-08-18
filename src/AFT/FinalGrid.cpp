
#include "AFT.h"

namespace AFT
{
    void addToPointADT (const Point& p, PointADT& pointADT, int index)
    {
        ADT::ADTPoint vec = pointADT.createADTPoint (p.dim, p.dim);
        vec.idx = index;
        bool tempBool;
        pointADT.insert (vec, pointADT.root, tempBool);
        if (tempBool == false)
        {
            cout << "could not insert to pointADT" << endl;
            exit(-2);
        }
    }
    
    bool cmpPoints (const ::Point& p0, const ::Point& p1)
    {
        int counter = 0;
        for (int d=0; d<N_DIM; ++d)
        {
            if ( fabs(p0.dim[d] - p1.dim[d]) < 1e-10 )
            {
                ++counter;
            }
        }

        if (counter == N_DIM)
        {
            return true;
        }
        
        return false;
    }
    
    bool pointExistsForCreateCells(const ::Point& refPoint, const vector<::Point>& points, int& index)
    {
        for (int ip=0; ip<points.size(); ++ip)
        {
            if ( cmpPoints (refPoint, points[ip]) )
            {
                index = ip;
                return true;
            }
        }
        
        return false;
    }
    
    bool cellExistsForCreateCells (const Cell& refCell, const vector<Cell>& cell, int& index)
    {
        for (int c=0; c<cell.size(); ++c)
        {
            int counter = 0;
            for (int d=0; d<N_DIM; ++d)
            {
                if ( fabs(refCell.cnt[d] - cell[c].cnt[d]) < 0.0001 )
                {
                    ++counter;
                }
            }

            if (counter == N_DIM)
            {
                index = c;
                return true;
            }
        }
        
        return false;
    }
    
    bool faceExistsForCreateCells (const Face& refFace, const vector<Face>& face, int& index)
    {
        for (int c=0; c<face.size(); ++c)
        {
            int counter = 0;
            for (int d=0; d<N_DIM; ++d)
            {
                if ( fabs(refFace.cnt[d] - face[c].cnt[d]) < 0.0001 )
                {
                    ++counter;
                }
            }
            
            /*if (refFace.cnt[0] < -6.98640 && refFace.cnt[0] > -6.98642)
            {
                if (face[c].cnt[0] < -6.98640 && face[c].cnt[0] > -6.98642)
                {
                    cout << refFace.cnt[0] << endl;
                    cout << refFace.cnt[1] << endl;
                    cout << refFace.cnt[2] << endl;
                    
                    cout << face[c].cnt[0] << endl;
                    cout << face[c].cnt[1] << endl;
                    cout << face[c].cnt[2] << endl;
                    
                    cout << counter << endl;
                    
                    for (int d=0; d<N_DIM; ++d)
                    {
                        if (refFace.cnt[d] == face[c].cnt[d])
                        {
                            cout << "equal_" << d << endl;
                        }
                        else
                        {
                            cout << refFace.cnt[d] << endl;
                            cout << face[c].cnt[d] << endl;
                        }
                    }
                    
                    cin.ignore();
                }
            }*/

            if (counter == N_DIM)
            {
                index = c;
                return true;
            }
        }
        
        return false;
    }
    
    bool strictlyInteriorPass (const Face& f, const iBlank_t& LCIBlank, const iBlank_t& RCIBlank, const bool LCTrim, const bool RCTrim)
    {
        if ( f.bouType != face_t::INTERIOR )
        {
            return false;
        }
        else if ( f.nei.size() != 2 )
        {
            return false;
        }
        else if (LCIBlank != iBlank_t::FIELD || RCIBlank != iBlank_t::FIELD)
        {
            return false;
        }
        else if (LCTrim != false || RCTrim != false)
        {
            return false;
        }        
        
        return true;
    }
    
    bool ghostPass (const Face& f, const iBlank_t& LCIBlank, const bool LCTrim)
    {
        if (f.bouType != face_t::BOUNDARY)
        {
            return false;
        }
        else if (LCIBlank != iBlank_t::FIELD || LCTrim != false)        
        {
            return false;
        }
        
        return true;
    }
    
    bool intergridPass (const Face& f, const iBlank_t& LCIBlank, const iBlank_t& RCIBlank,
            const bool LCTrim, const bool RCTrim)
    {
        if (f.bouType != face_t::INTERIOR)
        {
            return false;
        }
        else if (LCIBlank != iBlank_t::FIELD || RCIBlank != iBlank_t::FIELD) // both sides must be field
        {
            return false;
        }
        /*else if (LCTrim == false && RCTrim == false)
        {
            return false;
        }
        else if (LCTrim == true && RCTrim == true)
        {
            return false;
        }*/
        else if (LCTrim == false)
        {
            if (RCTrim == false)
            {
                return false;
            }
        }
        else if (LCTrim == true)
        {
            if (RCTrim == true)
            {
                return false;
            }
        }
                
        return true;
    }
    
    void addNearBoundaryCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<::Point>& pt, Grid& finalGrid, PointADT& fgp, PointADT& fgcc, PointADT& fgfc)
    {
        // face f is the face of normal or new grid
    
        int iRefC;
        int iRefFace;
        
        Point cellCent;
        cellCent.dim = cell[f.nei[0]].cnt;
        bool cellExist = pointExists (cellCent, fgcc, iRefC); // cell should exist
        if (!cellExist)
        {
            cout << "cell opps" << endl;
            exit(-2);
        }
        Cell& cll = finalGrid.cell[iRefC];
        
        Point faceCent;
        faceCent.dim = f.cnt;
        bool faceExist = pointExists (faceCent, fgfc, iRefFace); // face should exist
        Face& fac = finalGrid.face[iRefFace];
        
        cll.face.push_back (iRefFace);
        fac.nei.push_back (iRefC);
        
        int tmp = fac.nei[0];
        fac.nei[0] = fac.nei[1];
        fac.nei[1] = tmp;
        
        finalGrid.cell[fac.nei[0]].nei.push_back (fac.nei[1]);
        finalGrid.cell[fac.nei[1]].nei.push_back (fac.nei[0]);

        // loop through vertices of cell
        for (int v=0; v<f.vtx.size(); ++v)
        {
            for (int cv=0; cv<cll.vtx.size(); ++cv)
            {
                if (cll.vtxBelo[cv] != finalGrid.id)
                {
                    if (cll.vtx[cv] == f.vtx[v])
                    {
                        cll.vtx[cv] = fac.vtx[v];
                        cll.vtxBelo[cv] = finalGrid.id;
                        break;
                    }
                }
            }
        }
    }
    
    void addGhostCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<::Point>& pt, Grid& finalGrid, PointADT& fgp, PointADT& fgcc, PointADT& fgfc)
    {
        // create new face
        Face newFace = f;
        newFace.vtx.clear();
        newFace.nei.clear();
        int iRefFace = finalGrid.face.size();
     
        // create new cell
        Cell newCell = cell[f.nei[1]]; // boundary element
        newCell.nei.clear();
        newCell.face.clear();
        //newCell.belonging = finalGrid.id;
        int iRefC = finalGrid.cell.size();
        Point cellCent;
        cellCent.dim = newCell.cnt;
        addToPointADT (cellCent, fgcc, iRefC);
        newCell.face.push_back(iRefFace);
        finalGrid.cell.push_back (move(newCell));
        Cell& cll = finalGrid.cell [iRefC];
        newFace.nei.push_back (iRefC);        

        // loop through vertices of face
        for (int v=0; v<f.vtx.size(); ++v)
        {
            // create new point
            ::Point p = pt[f.vtx[v]];
            Point pp;
            pp.dim = p.dim;
            
            int pointIndex;
            bool ptExists = pointExists (pp, fgp, pointIndex);

            if (!ptExists)
            {
                pp.belonging = finalGrid.id;
                finalGrid.pt.push_back (move(p));
                pointIndex = finalGrid.pt.size() - 1;
                addToPointADT (pp, fgp, pointIndex);
            }

            newFace.vtx.push_back( pointIndex );
            
            // loop through vertices of cell
            for (int cv=0; cv<cll.vtx.size(); ++cv)
            {
                if (cll.vtxBelo[cv] != finalGrid.id)
                {
                    if (cll.vtx[cv] == f.vtx[v])
                    {
                        cll.vtx[cv] = pointIndex;
                        cll.vtxBelo[cv] = finalGrid.id;
                        break;
                    }
                }
            }
        }
        
        Point faceCent;
        faceCent.dim = newFace.cnt;
        addToPointADT (faceCent, fgfc, iRefFace);
        finalGrid.face.push_back (move(newFace));
    }
        
    void addCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<::Point>& pt, Grid& finalGrid, PointADT& fgp, PointADT& fgcc)
    {
        int neiSize = f.nei.size();
        if (neiSize != 2)
        {
            cout << "neiSize != 2" << endl;
            exit(-2);
        }
        int iRefC;
    
        // create new face
        Face newFace = f;
        newFace.vtx.clear();
        newFace.nei.clear();
        int iRefFace = finalGrid.face.size();
        
        for (int i=0; i<neiSize; ++i)
        {
            Point cellCent;
            cellCent.dim = cell[f.nei[i]].cnt;
            bool cellExist = pointExists (cellCent, fgcc, iRefC);
            
            if (!cellExist)
            {
                Cell newCell = cell[f.nei[i]];
                newCell.nei.clear();
                newCell.face.clear();
                //newCell.belonging = finalGrid.id;
                finalGrid.cell.push_back (move(newCell));
                iRefC = finalGrid.cell.size() - 1;
                addToPointADT (cellCent, fgcc, iRefC);
            }
            
            finalGrid.cell[iRefC].face.push_back(iRefFace);
            newFace.nei.push_back (iRefC);
        }
        
        if (newFace.nei.size() == 2)
        {
            finalGrid.cell[newFace.nei[0]].nei.push_back (newFace.nei[1]);
            finalGrid.cell[newFace.nei[1]].nei.push_back (newFace.nei[0]);
        }
        else
        {
            cout << "should be two in finalgrid" << endl;
            exit(-2);
        }
        
        // loop through vertices of face
        for (int v=0; v<f.vtx.size(); ++v)
        {
            // create new point
            ::Point p = pt[f.vtx[v]];
            Point pp;
            pp.dim = p.dim;

            int pointIndex = -1;
            bool ptExists = pointExists (pp, fgp, pointIndex);

            if (!ptExists)
            {
                pp.belonging = finalGrid.id;
                finalGrid.pt.push_back (p);
                pointIndex = finalGrid.pt.size() - 1;
                addToPointADT (pp, fgp, pointIndex);
            }

            newFace.vtx.push_back( pointIndex );
            
            // loop through vertices of cell
            for (int i=0; i<neiSize; ++i)
            {
                Cell& cll = finalGrid.cell[newFace.nei[i]];
            
                for (int cv=0; cv<cll.vtx.size(); ++cv)
                {
                    if (cll.vtxBelo[cv] != finalGrid.id)
                    {
                        if (cll.vtx[cv] == f.vtx[v])
                        {
                            cll.vtx[cv] = pointIndex;
                            cll.vtxBelo[cv] = finalGrid.id;
                            break;
                        }
                    }
                }
            }
        }
        
        finalGrid.face.push_back (move(newFace));
    }
    
    void modifyCellVertices (Grid& finalGrid, const Grid& newGrid, const vector<Grid>& gr, PointADT& fgp)
    {
        bool match;
        int pointIndex;
        
        for (int ic=finalGrid.n_bou_elm; ic<finalGrid.cell.size(); ++ic)
        {
            Cell& c = finalGrid.cell[ic];
            
            for (int& v: c.vtx)
            {
                if (c.belonging == 0 || c.belonging == 1)
                {
                    Point pp;
                    pp.dim = gr[c.belonging].pt[v].dim;
                    match = pointExists (pp, fgp, pointIndex);
                }
                else if (c.belonging == 2)
                {
                    Point pp;
                    pp.dim = newGrid.pt[v].dim;
                    match = pointExists (pp, fgp, pointIndex);
                }
                else
                {
                    cout << "!!! Error: c.belonging != 0 || c.belonging != 1 || c.belonging != 2 in modifyCellVertices(...)" << endl;
                    cout << "!!! Error: c.belonging = " << c.belonging << endl;
                    exit(-2);
                }
                
                if (match)
                {
                    v = pointIndex;
                }
                else
                {
                    cout << "!!! Error: no match in modifyCellVertices(...)" << endl;
                    cout << "!!! Error: c.belonging = " << c.belonging << endl;
                    cout << "!!! Error: dim[0] = " << gr[c.belonging].pt[v].dim[0] << endl;
                    cout << "!!! Error: dim[1] = " << gr[c.belonging].pt[v].dim[1] << endl;
                    cout << "!!! Error: dim[2] = " << gr[c.belonging].pt[v].dim[2] << endl;
                    exit(-2);
                }
            }
            
            c.belonging = finalGrid.id;
        }
    }
    
    void addIntergridCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<::Point>& pt, Grid& finalGrid, const Grid& newGrid, PointADT& fgp, PointADT& fgcc, PointADT& ngfc)
    {
        bool accorL = true;
        bool cellExist;
    
        int neiSize = f.nei.size();
        if (neiSize != 2)
        {
            cout << "neiSize != 2" << endl;
            exit(-2);
        }
        int iRefC;
    
        // create new face
        Face newFace = f;
        newFace.vtx.clear();
        newFace.nei.clear();
        int iRefFace = finalGrid.face.size();
        
        // find and add existing normal grid cell
        Point cellCent;
        if (cell[f.nei[0]].trim == false)
        {
            cellCent.dim = cell[f.nei[0]].cnt;
            cellExist = pointExists (cellCent, fgcc, iRefC);
        }
        else
        {
            // face conf should be changed for right cell
            accorL = false;
            cellCent.dim = cell[f.nei[1]].cnt;
            cellExist = pointExists (cellCent, fgcc, iRefC);
        }
        
        if (!cellExist)
        {
            cout << "cell must exist 1 in add intergrid" << endl;
            exit(-2);
        }
        
        newFace.nei.push_back (iRefC);
        finalGrid.cell[iRefC].face.push_back (iRefFace);
        
        // find and add existing new grid cell
        int faceIndex;
        Point faceCent;
        faceCent.dim = newFace.cnt;
        pointExists (faceCent, ngfc, faceIndex);
        const Face& ngf = newGrid.face[faceIndex];
        cellCent.dim = newGrid.cell [ngf.nei[0]].cnt;
        cellExist = pointExists (cellCent, fgcc, iRefC); // should exist
        if (!cellExist)
        {
            cout << "cell must exist in add intergrid" << endl;
            exit(-2);
        }
        
        finalGrid.cell [iRefC].face.push_back (iRefFace);
        newFace.nei.push_back (iRefC);

        finalGrid.cell[newFace.nei[0]].nei.push_back(newFace.nei[1]);
        finalGrid.cell[newFace.nei[1]].nei.push_back(newFace.nei[0]);

        // loop through vertices of normal grid face
        for (int v=0; v<f.vtx.size(); ++v)
        {
            // create new point
            ::Point p = pt[ f.vtx[v] ];
            Point pp;
            pp.dim = p.dim;

            int pointIndex;
            bool ptExists = pointExists (pp, fgp, pointIndex);

            if (!ptExists)
            {
                pp.belonging = finalGrid.id;
                finalGrid.pt.push_back (p);
                pointIndex = finalGrid.pt.size() - 1;
                addToPointADT (pp, fgp, pointIndex);
            }

            newFace.vtx.push_back( pointIndex );
            
            // loop through vertices of cell            
            Cell& cll = finalGrid.cell[newFace.nei[0]];
        
            for (int cv=0; cv<cll.vtx.size(); ++cv)
            {
                if (cll.vtxBelo[cv] != finalGrid.id)
                {
                    if (cll.vtx[cv] == f.vtx[v])
                    {
                        cll.vtx[cv] = pointIndex;
                        cll.vtxBelo[cv] = finalGrid.id;
                        break;
                    }
                }
            }
        }
        
        // loop through vertices of new grid face
        for (int v=0; v<ngf.vtx.size(); ++v)
        {
            // create new point
            ::Point p = newGrid.pt[ ngf.vtx[v] ];
            Point pp;
            pp.dim = p.dim;

            int pointIndex;
            bool ptExists = pointExists (pp, fgp, pointIndex); // should exist
        
            // loop through vertices of cell            
            Cell& cll = finalGrid.cell[newFace.nei[1]];
        
            for (int cv=0; cv<cll.vtx.size(); ++cv)
            {
                if (cll.vtxBelo[cv] != finalGrid.id)
                {
                    if (cll.vtx[cv] == ngf.vtx[v])
                    {
                        cll.vtx[cv] = pointIndex;
                        cll.vtxBelo[cv] = finalGrid.id;
                        break;
                    }
                }
            }
        }
        
        if (!accorL)
        {
            int tmp = newFace.vtx[1];
            newFace.vtx[1] = newFace.vtx.back();
            newFace.vtx.back() = tmp;
            newFace.set_area (finalGrid.pt);
        }
        finalGrid.face.push_back (move(newFace));
    }
    
    void findOtherNeiOfGhostsFaces (Grid& finalGrid, PointADT& fgcc)
    {
        for (Face& f: finalGrid.face)
        {
            if (f.bouType == face_t::BOUNDARY)
            {
                if (f.nei.size() == 1)
                {
                    Point cellCent;
                    cellCent.dim = f.cnt;
                    int ic;                    
                    bool pExists = pointExists (cellCent, fgcc, ic);
                    
                    if (pExists)
                    {
                        Cell& cll = finalGrid.cell[ic];
                        
                        f.nei.push_back(ic);
                        finalGrid.cell[f.nei[0]].nei.push_back(ic);
                        cll.nei.push_back(f.nei[0]);
                        int tmp = f.nei[0];
                        f.nei[0] = f.nei[1];
                        f.nei[1] = tmp;

                        cll.face.push_back (&f - &finalGrid.face.front());
                    }
                    else
                    {
                        cout << "could not find in findOtherNeiOfGhostsFaces (...)" << endl;
                        exit(-2);
                    }
                }
                else
                {
                    cout << "!!! Error: f.nei.size() != 1 in findOtherNeiOfGhostsFaces(...)" << endl;
                    cout << "!!! Error: f.nei.size() = " << f.nei.size() << endl;
                    exit(-2);
                }
            }
        }
    }
    
    void buildPointADTforFinalGrid (PointADT& pADT, const vector<Grid>& gr)
    {
        pADT.root = new PointADT::Node();
        pADT.root->level = 0;
        //pADT.root->p->idx = -1;
        
        for (int d=0; d<ADT_VAR; ++d)
        {
            pADT.root->c[d] = BIG_POS_NUM;
            pADT.root->d[d] = BIG_NEG_NUM;
        }
        
        for (const Grid& g: gr)
        {
            for (const ::Point& p: g.pt)
            {
                for (int d=0; d<ADT_DIM; ++d)
                {
                    pADT.root->c[2*d] = min (p.dim[d], pADT.root->c[2*d]);
                    pADT.root->d[2*d] = max (p.dim[d], pADT.root->d[2*d]);
                }
            }
        }
        
        for (int d=0; d<ADT_DIM; ++d)
        {
            pADT.root->c[2*d+1] = pADT.root->c[2*d];
            pADT.root->d[2*d+1] = pADT.root->d[2*d];
        }
    }
    
    void createFinalGrid (Grid& finalGrid, const vector<Grid>& gr, const Grid& newGrid)
    {
        // reserve finalgrid.cell!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        PointADT fgp; // final grid points
        PointADT fgcc; // final grid cell centers
        PointADT fgfc; // final grid boundary face centers
        PointADT ngfc; // new grid face centers
        
        //cout << "Building fgp..." << flush;
        buildPointADTforFinalGrid (fgp, gr);
        //cout << " done" << endl;
        
        //cout << "Building fgcc..." << flush;
        buildPointADTforFinalGrid (fgcc, gr);
        //cout << " done" << endl;
        
        //cout << "Building fgcc..." << flush;
        buildPointADTforFinalGrid (fgfc, gr);
        //cout << " done" << endl;
        
        
        //cout << "Building ngfc..." << flush;
        ngfc.points.reserve (newGrid.face.size());
        
        for (int iF=0; iF<newGrid.face.size(); ++iF)                
        {
            const Face& f = newGrid.face[iF];
            
            ngfc.points.push_back ( ngfc.createADTPoint (f.cnt, f.cnt) );
            ngfc.points.back().idx = iF;
        }
        
        ngfc.build();
        
        
        cout << "test in final grid" << endl;
        
        // add boundary elements which are attached to field and untrimmed cells        
        for (const Grid& g: gr)
        {
            for (int ff=0; ff<g.face.size(); ++ff)
            //for (const Face& f: g.face)
            {
                const Face& f = g.face[ff];
            
                bool pass = ghostPass (f, g.cell[f.nei[0]].iBlank, g.cell[f.nei[0]].trim);
                if (pass) {addGhostCells (f, g.face, g.cell, g.pt, finalGrid, fgp, fgcc, fgfc);}
            }
        }
        
        for (const Face& f: newGrid.face)
        {
            bool pass = ghostPass (f, newGrid.cell[f.nei[0]].iBlank, newGrid.cell[f.nei[0]].trim);
            if (pass) { addGhostCells (f, newGrid.face, newGrid.cell, newGrid.pt, finalGrid, fgp, fgcc, fgfc);}
        }
        
        finalGrid.n_bou_elm = finalGrid.cell.size();
        
        // add strictly interior cells of normal grids        
        for (int gg=0; gg<gr.size(); ++gg)
        {
            const Grid& g = gr[gg];
        
            for (int ff=0; ff<g.face.size(); ++ff)            
            {
                const Face& f = g.face[ff];
                
                bool pass = strictlyInteriorPass (f, g.cell[f.nei[0]].iBlank, g.cell[f.nei[1]].iBlank,
                                                  g.cell[f.nei[0]].trim, g.cell[f.nei[1]].trim);
                if (pass)
                {
                    addCells (f, g.face, g.cell, g.pt, finalGrid, fgp, fgcc);
                }
            }
        }
        
        //cout << "finalGrid.cell.size() = " << finalGrid.cell.size() << endl;
        
        // add strictly interior cells of new grid
        for (const Face& f: newGrid.face)
        {
            if (f.nei.size() == 2)
            {
                bool pass = strictlyInteriorPass (f, newGrid.cell[f.nei[0]].iBlank, newGrid.cell[f.nei[1]].iBlank,
                                                  newGrid.cell[f.nei[0]].trim, newGrid.cell[f.nei[1]].trim);
                
                if (pass) addCells (f, newGrid.face, newGrid.cell, newGrid.pt, finalGrid, fgp, fgcc);
            }
        }
        
        //cout << "finalGrid.cell.size() = " << finalGrid.cell.size() << endl;
        
        // add intergrid cells        
        for (int gg=0; gg<gr.size(); ++gg)
        {
            const Grid& g = gr[gg];
        
            for (int ff=0; ff<g.face.size(); ++ff)
            {
                const Face& f = g.face[ff];
                
                bool pass = intergridPass (f, g.cell[f.nei[0]].iBlank, g.cell[f.nei[1]].iBlank,
                                                  g.cell[f.nei[0]].trim, g.cell[f.nei[1]].trim);
                
                if (pass)
                {
                    addIntergridCells (f, g.face, g.cell, g.pt, finalGrid, newGrid, fgp, fgcc, ngfc);
                }
            }
        }
        
        //cout << "finalGrid.cell.size() = " << finalGrid.cell.size() << endl;
        
        // connect ghosts and interiors
        for (const Grid& g: gr)
        {
            for (const Face& f: g.face)
            {
                bool pass = ghostPass (f, g.cell[f.nei[0]].iBlank, g.cell[f.nei[0]].trim);
                if (pass) {addNearBoundaryCells (f, g.face, g.cell, g.pt, finalGrid, fgp, fgcc, fgfc);}
            }
        }
        
        for (const Face& f: newGrid.face)
        {
            bool pass = ghostPass (f, newGrid.cell[f.nei[0]].iBlank, newGrid.cell[f.nei[0]].trim);
            if (pass) { addNearBoundaryCells (f, newGrid.face, newGrid.cell, newGrid.pt, finalGrid, fgp, fgcc, fgfc); }
        }
        
        finalGrid.totalNElms = finalGrid.cell.size();
        finalGrid.n_in_elm = finalGrid.totalNElms - finalGrid.n_bou_elm;
        
        //cout << "finalGrid.n_bou_elm = " << finalGrid.n_bou_elm << endl;
        //cout << "finalGrid.n_in_elm = " << finalGrid.n_in_elm << endl;
        //cout << "finalGrid.totalNElms = " << finalGrid.totalNElms << endl;
    }
}
