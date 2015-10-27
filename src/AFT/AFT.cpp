#include "AFT.h"

namespace AFT
{
    void aft (vector<Grid>& gr, Grid& finalGrid)
    {
        //cout << "Allocating... " << flush;
        int newGridId = 2;        
        double offsetZ = 1.;
        int phys = 1;
        Grid newGrid (gr[0].mainDir, newGridId);        
        double aveTriArea;
        vector<Point> points;
        vector<Point> edgeCenters;
        vector<Edge> edges;
        vector<Edge> mesh0Edges;
        vector<Edge> mesh1Edges;
        vector<Triangle> triangles;
        vector<FrontMember> frontList;
        EdgeADT edgeADT;
        EdgeADT edge0ADT;
        EdgeADT edge1ADT;
        EdgeADT edge01ADT;
        TriangleADT triangleADT;
        PointADT pointADT;
        PointADT edgeCenterADT;
        CircleADT circleADT;
        //cout << "done!" << endl;
        
        // this should not be in AFT
        //cout << "Trimming/Re-blanking... " << flush;
        gr[0].trimWhoHasFringeNeighbor();
        gr[1].trimWhoHasFringeNeighbor();
        gr[0].trimWhoHasTrimNeighbor (3);
        gr[1].trimWhoHasTrimNeighbor (2);
        //cout << "done!" << endl;
        
        // this too
        //cout << "Outputing after trimming/reblanking... " << flush;
        //gr[0].outAllTecplot();
        //gr[1].outAllTecplot();
        gr[0].outAllVTK(0);
        gr[1].outAllVTK(0);
        //cout << "done!" << endl;
     
        //cout << "Preparing... " << flush;
        setPointsEdges (gr, points, edges, edgeCenters, newGridId);
        createFrontList (edges, frontList, points);
        aveTriArea = getAveTriArea (edges, points);
        //cout << "done!" << endl;
        
        //cout << "Building trees... " << flush;        
        for (Edge& e: edges)
        {
            if (e.belonging == 0)
            {
                mesh0Edges.push_back (e);
            }
            else if (e.belonging == 1)
            {
                mesh1Edges.push_back (e);
            }
        }
        
        edge0ADT.build (points, mesh0Edges);
        edge1ADT.build (points, mesh1Edges);
        edgeADT.build (points, edges);
        edge01ADT.build (points, edges);
        triangleADT.build (edgeADT);
        pointADT.build (points);
        edgeCenterADT.build (edgeCenters);
        circleADT.build (edgeADT);
        //cout << "done!" << endl;
        
        //cout << "Exporting to GMSH... " << flush;
        exportToGMSH (points, mesh0Edges, mesh1Edges, gr[0].mainDir);
        //cout << "done!" << endl;        
        
        //cout << "Advancing front... " << flush;
        advanceFront (frontList, points, aveTriArea, edges, triangles, triangleADT, pointADT, edgeCenterADT, edgeADT, edge01ADT, newGridId, edgeCenters, circleADT);
        //cout << "done!" << endl;
        
        //cout << "Erasing dead elements... " << flush;
        eraseDeadPoints (points, edges, triangles);
        eraseDeadEdges (edges, triangles, points);
        eraseDeadTriangles (triangles, points, edges);
        //cout << "done!" << endl;
        
        //cout << "Outputing unflipped triangles... " << flush;        
        outputTrianglesVTK (points, triangles, gr[0].mainDir, "tri.vtk");
        //cout << "done!" << endl;        
        
        //cout << "Flipping triangles... " << flush;
        flip (triangles, edges, points);
        //cout << "done!" << endl;
        
        //cout << "Flipping triangles... " << flush;
        flip (triangles, edges, points);
        //cout << "done!" << endl;
        
        //cout << "Outputing flipped triangles... " << flush;    
        outputTrianglesVTK (points, triangles, gr[0].mainDir, "triFlip.vtk");
        //cout << "done!" << endl;
        
        //cout << "Creating cells... " << flush;
        createCells (offsetZ, points, newGrid, triangles, phys, newGridId);
        //cout << "done!" << endl;
        
        
        
        
        //cout << "Outputing new grid... " << flush;
        //newGrid.outAllTecplot();
        newGrid.outAllVTK(0);
        //cout << "done!" << endl;
        
        //cout << "Creating final grid... " << flush;
        createFinalGrid (finalGrid, gr, newGrid);
        //cout << "done!" << endl;
        
        
        
        //cout << "Outputing final grid... " << flush;
        //finalGrid.outAllTecplot();
        //finalGrid.outAllVTK(0);
        //cout << "done!" << endl;
        
        cout << "finished AFT" << endl;
    }
    
    void construct (int iCPX, bool A_CPX_exists, bool B_CPX_exists, int iA_CPX, int iB_CPX, int iA, int iB, vector<FrontMember>& frontList,
             vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, vector<Point>& points, CircleADT& circleADT, int iFrontEdge)
    {
        //int iFrontEdge = frontList.front().edge;
        
        if (!A_CPX_exists)
        {
            Edge tmpEdge = createEdge (iA, iCPX, newGridId, true);
            addToEdgeList (tmpEdge, iA, iCPX, edges, edgeADT, points);
            iA_CPX = edges.size() - 1;
            addToFrontList (iA_CPX, frontList, points, edges);
            
            Point cntPoint;
            cntPoint.belonging = newGridId;
            cntPoint.dim = 0.5 * (points[tmpEdge.t[0]].dim + points[tmpEdge.t[1]].dim);
            //addToPointList (cntPoint, edgeCenters, edgeCenterADT);
        }
        else
        {
            eraseExistingEdgeFromFrontList (iA_CPX, frontList);
        }
        
        if (!B_CPX_exists)
        {
            Edge tmpEdge = createEdge (iB, iCPX, newGridId, true);
            addToEdgeList (tmpEdge, iB, iCPX, edges, edgeADT, points);
            iB_CPX = edges.size() - 1;
            addToFrontList (iB_CPX, frontList, points, edges);
            
            Point cntPoint;
            cntPoint.belonging = newGridId;
            cntPoint.dim = 0.5 * (points[tmpEdge.t[0]].dim + points[tmpEdge.t[1]].dim);
            //addToPointList (cntPoint, edgeCenters, edgeCenterADT);
        }
        else
        {
            eraseExistingEdgeFromFrontList (iB_CPX, frontList);
        }
        
        Triangle tmpTriangle = createTriangle (iFrontEdge, iA_CPX, iB_CPX, edges, points);
        addToTriangleList (triangles, tmpTriangle, triangleADT, points, circleADT, edges, frontList);
        
        eraseFromFrontList (frontList, iFrontEdge);
        //sortFrontList (frontList, points, edges);
    }
    
    void exportToGMSH (const vector<Point>& points, const vector<Edge>& mesh0Edges, const vector<Edge>& mesh1Edges, string dir)
    {
        ofstream out;
        dir.append ("/exprtGMSH.geo");
        out.open (dir.c_str());
        
        vector <vector <int> > ptConn ((points.size()+1));
        
        if (out.is_open())
        {
            out << "Mesh.Algorithm = 6;" << endl;
            
            for (int i=0; i<points.size(); ++i)
            {
                out << "Point(";
                out << i;
                out << ") = {";
                out << points[i].dim[0];
                out << ",";
                out << points[i].dim[1];
                out << ",";
                out << points[i].dim[2];
                out << "};" << endl;
            }
            
            //---------------------------------------------
            {
                int start, end;
                
                start = 0;
                end = start + mesh0Edges.size();

                for (int e=start; e<end; ++e) // e for edge
                {
                    int p0 = mesh0Edges[e-start].t[0] ;
                    int p1 = mesh0Edges[e-start].t[1];

                    ptConn[p0].push_back(e);
                    ptConn[p1].push_back(e);
                }

                start = mesh0Edges.size();
                end = start + mesh1Edges.size();

                for (int e=start; e<end; ++e)
                {
                    int p0 = mesh1Edges[e-start].t[0];
                    int p1 = mesh1Edges[e-start].t[1];

                    ptConn[p0].push_back(e);
                    ptConn[p1].push_back(e);
                }
            }
            //---------------------------------------------
            
            out << "Line(";
            out << 0;
            out << ") = {";
            out << mesh0Edges[0].t[0];            
            out << ",";
            out << mesh0Edges[0].t[1];
            out << "};" << endl;
            
            int last = mesh0Edges[0].t[1];
            int iLastEdge = 0;
            
            for (int i=0; i<mesh0Edges.size()-1; ++i)
            {
                // find the edge that last belongs
                for (int j=0; j<2; ++j)
                {
                    if ( ptConn[last][j] != iLastEdge )
                    {
                        // found j for iOtherEdge
                        int iOtherEdge = ptConn[last][j];
                        iLastEdge = iOtherEdge;             
                        
                        // find iOtherEdge.t != it1
                        for (int k=0; k<2; ++k)
                        {
                            if ( mesh0Edges[iOtherEdge].t[k] != last )
                            {                                
                                // found k for iOtherEdge.t != it1                                
                                out << "Line(";
                                out << i+1;
                                out << ") = {";
                                out << last;
                                out << ",";
                                out << mesh0Edges[iOtherEdge].t[k];
                                out << "};" << endl;
                                
                                last = mesh0Edges[iOtherEdge].t[k];                                
                                break;
                            }
                        }
                        
                        break;
                    }
                }
            }
            
            out << "Line(";
            out << mesh0Edges.size();
            out << ") = {";
            out << mesh1Edges[0].t[0];
            out << ",";
            out << mesh1Edges[0].t[1];
            out << "};" << endl;
            
            last = mesh1Edges[0].t[1];
            iLastEdge = mesh0Edges.size();
            
            for (int i=0; i<mesh1Edges.size()-1; ++i)
            {
                // find the edge that last belongs
                for (int j=0; j<2; ++j)
                {
                    if ( ptConn[last][j] != iLastEdge )
                    {
                        // found j for iOtherEdge
                        int iOtherEdge = ptConn[last][j] - mesh0Edges.size();
                        iLastEdge = ptConn[last][j];             
                        
                        // find iOtherEdge.t != it1
                        for (int k=0; k<2; ++k)
                        {
                            if ( mesh1Edges[iOtherEdge].t[k] != last )
                            {
                                // found k for iOtherEdge.t != it1                                
                                out << "Line(";
                                out << i+1 + mesh0Edges.size();
                                out << ") = {";
                                out << last;
                                out << ",";
                                out << mesh1Edges[iOtherEdge].t[k];
                                out << "};" << endl;
                                
                                last = mesh1Edges[iOtherEdge].t[k];
                                break;
                            }
                        }
                        
                        break;
                    }
                }
            }
            
            out << "Line Loop(1) = {";
            for (int i=0; i<mesh0Edges.size(); ++i)
            {
                out << i;
                
                if (i < mesh0Edges.size()-1)
                {
                    out << ",";
                }
            }
            out << "};" << endl;
            
            out << "Line Loop(2) = {";
            for (int i=mesh0Edges.size(); i<mesh0Edges.size()+mesh1Edges.size(); ++i)
            {
                out << i;
                
                if (i < mesh0Edges.size()+mesh1Edges.size()-1)
                {
                    out << ",";
                }
            }
            out << "};" << endl;
            
            out << "Plane Surface(1) = {1,2};";
            
            out.close();
        }
        else
        {
            cout << "could not open file in AFT::exportToGMSH()" << endl;
            exit(-2);
        }
    }
    
    double spacingFnc (double b, double aveTriSize)
    {
        // initially assume uniform mesh
        
        double h = 2. * aveTriSize / b;
        
        return (h/2.);
    }
    
    void Tanemura_Merriam_Helper (int iA, int iB, int& i, vector<Point>& points, deque<int>& pts)
    {
        bool found = false;
        
        CVector center;
        double radius;
        
        Point& A = points[iA];
        Point& B = points[iB];
        
        Point& CPX = points[pts[i]];
        triPtsCircums (CPX.dim, A.dim, B.dim, center, radius);        
                
        //cout << "iA = " << iA << endl;
        //cout << "iB = " << iB << endl;
        //cout << "pts.size() = " << pts.size() << endl;
        
        for (int j=0; j<pts.size(); ++j)
        {
            //cout << "j = " << j << endl;
            //cout << "i = " << i << endl;
            //cout << "pts[j] = " << pts[j] << endl;
            if (j != i && pts[j] != -1)
            {
                CVector d = points[pts[j]].dim - center;

                if ( (radius - mag(d)) > 1e-5 )
                {
                    pts[i] = -1;
                    i = j;
                    //cout << "pts[i] = " << pts[i] << endl;                    
                    //pts.erase (pts.begin() + j);
                    //cout << "radius = " << radius << endl;
                    //cout << "mag(d) = " << mag(d) << endl;
                    found = true;
                    break;
                }
                else
                {
                    pts[j] = -1;
                }
            }
        }
                    
        if ( found )
        {
            Tanemura_Merriam_Helper (iA, iB, i, points, pts);
        }
    }
    
    int Tanemura_Merriam (int iA, int iB, vector<Point>& points, deque<int> pts)
    {
        //cout << "iA = " << iA << endl;
        //cout << "iB = " << iB << endl;        
    
        if (pts.size() == 0)
        {
            cout << "pts.size() == 0 in AFT::Tanemura_Merriam(...)" << endl;
            return -1;
            //exit(-2);
        }
        else
        {
            int i = -1;
        
            for (int j=0; j<pts.size(); ++j)
            {
                //cout << "pts[j] = " << pts[j] << endl;
            
                if (pts[j] == iA || pts[j] == iB)
                {
                    pts[j] = -1;
                }
            }
            
            for (int j=0; j<pts.size(); ++j)
            {
                if (pts[j] != -1)
                {
                    i = j;
                    break;
                }
            }
            
            if (i == -1)
            {
                cout << "no good points in AFT::Tanemura_Merriam(...)" << endl;
                return -1;
            }
        
            Tanemura_Merriam_Helper (iA, iB, i, points, pts);
        
            return i;
        }
    }
}
