//#include <bits/stl_vector.h>

#include "AFT.h"

namespace AFT
{
    FrontMember::FrontMember()
    {
        newPointChecked = false;
        //cloPtsMaxSize = 30;
    }
    
    void sortFrontList (vector<FrontMember>& frontList, const vector<Point>& points, const vector<Edge>& edges)
    {
        auto cmp = [&] (const FrontMember& fm1, const FrontMember& fm2)
        {
            //cout << "fm2.edge = " << fm2.edge << endl;
            
            const Edge& e1 = edges[ fm1.edge ];
            const Edge& e2 = edges[ fm2.edge ];
            
            const Point& e1t0 = points[ e1.t[0] ];
            const Point& e1t1 = points[ e1.t[1] ];
            const Point& e2t0 = points[ e2.t[0] ];
            const Point& e2t1 = points[ e2.t[1] ];
            
            double mag10 = pow( e1t0.dim[0] - e1t1.dim[0], 2 );
            double mag11 = pow( e1t0.dim[1] - e1t1.dim[1], 2 );
            double mag1  = mag10 + mag11;

            double mag20 = pow( e2t0.dim[0] - e2t1.dim[0], 2 );
            double mag21 = pow( e2t0.dim[1] - e2t1.dim[1], 2 );
            double mag2  = mag20 + mag21;

            if (e1.belonging < e2.belonging)
            {
                return true;
            }
            else if (e1.belonging > e2.belonging)
            {
                return false;
            }
            else if (mag1 < mag2)
            {
                return true;
            }

            return false;
        };
        
        sort(frontList.begin(), frontList.end(), cmp);
        
        /*for (int i=0; i<frontList.size(); ++i)
        {
            cout << edges[frontList[i].edge].belonging << endl;
        }
        
        cin.ignore();*/
    }
    
    void eraseFromFrontList (vector<FrontMember>& frontList, int edge)
    {
        if (frontList.size() != 0)
        {
            frontList.erase (remove_if (frontList.begin(), frontList.end(), [&](FrontMember& fm) { return fm.edge == edge; }), frontList.end());
        }
    }
    
    void addToFrontList (int edge, vector<FrontMember>& frontList, vector<Point>& points, vector<Edge>& edges)
    {
        FrontMember fm;
        fm.edge = edge;
        
        frontList.push_back(fm);
        
        sortFrontList (frontList, points, edges);
    }
    
    void eraseExistingEdgeFromFrontList (int ie, vector<FrontMember>& frontList)
    {
        for (unsigned int i=0; i<frontList.size(); ++i)
        {
            if ( frontList[i].edge == ie )
            {
                frontList.erase (frontList.begin() + i);
                return;
            }
        }
    }

    bool constructWithExisting (bool& foundCandPts, int& iFrontEdge, bool& existingPointPass, bool& A_CPX_exists, bool& B_CPX_exists, int& iA_CPX, int& iB_CPX, bool& edgesPass, const Point& t0, const Point& t1, bool& circumBoundPass, Point& CPX, int& iCPX, int& it0, int& it1, double& rho, deque<int>& candPts, FrontMember& frontFirst, vector<FrontMember>& frontList, vector<Point>& points, 
                       vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT,
                       PointADT& pointADT, EdgeADT& edgeADT, EdgeADT& edge01ADT, int newGridId, CircleADT& circleADT, int countr)
    {
		// search existing candidate points
		foundCandPts = srchCandPts (frontFirst, edges, points, pointADT, candPts, (2.*rho), edgeADT, edge01ADT, triangleADT);

		if (candPts.size() == 0)
		{
			outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
			//existingPointPass = false;
			//exit(-2);
		}

		if (foundCandPts)
		{
			// find Delaunay triangle
			outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
			int iCandPt = Tanemura_Merriam (it0, it1, points, candPts);
			iCPX = candPts[iCandPt];
			CPX = points[iCPX];

			// check whether forming triangle has a circumradius smaller than threshold radius
            if (CPX.dim[0] == t0.dim[0])
            {
                assert(CPX.dim[1] != t0.dim[1]);
            }
            if (CPX.dim[0] == t1.dim[0])
            {
                assert(CPX.dim[1] != t1.dim[1]);
            }
			circumBoundPass = checkCircumBound (CPX, t0, t1, rho);

			// if circumradius check passed check two forming edges intersections
			if (circumBoundPass)
			{
				edgesPass = checkTwoFormingEdges (CPX, t0, t1, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, edges, edgeADT, points);
			}

			existingPointPass = edgesPass;

            //if (iFrontEdge == 502)
            //{
            //    for (auto aa: candPts)
            //    {
            //        cout << "cp: " << aa << endl;
            //    }
            //    cout << "iCPX: " << iCPX << endl;
            //    outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
            //    assert(false);
            //}

            //if (countr == 36)
            {
                //if (it0 == 166 && it1 == 208)
                //if (it0 == 208 && it1 == 211)
                //{
                //    for (auto aa: candPts)
                //    {
                //        cout << "cp: " << aa << endl;
                //    }
                //    cout << "iCPX: " << iCPX << endl;
                //    cout << "circumBoundPass: " << circumBoundPass << endl;
                //    cout << "existingPointPass: " << existingPointPass << endl;
                //    outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
                //    assert(false);
                //}
            }
		}

		//cout << "iCPX = " << iCPX << endl;

		// construct triangle if all previous stages are passed
		if (existingPointPass)
		{
			//cout << "with exist point" << endl;

			construct (iCPX, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, it0,
					it1, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points, circleADT, iFrontEdge);

                    //for (int ppp = 0; ppp < points.size(); ++ppp)
                    //{
                    //    auto& pp = points[ppp];
                    //    if ((pp.dim[0] == CPX.dim[0] && pp.dim[1] == CPX.dim[1]))
                    //    {
                    //        if (ppp != iCPX)
                    //        {
                    //            assert(false);
                    //        }
                    //    }
                    //}



			/*if (triangles.size() == 164)
			  {

			  cout << "it0 = " << it0 << endl;
			  cout << "it1 = " << it1 << endl;
			  cout << "iCPX = " << iCPX << endl;
			  cout << "A_CPX_exists = " << A_CPX_exists << endl;
			  cout << "B_CPX_exists = " << B_CPX_exists << endl;
			  exit(-2);
			  }*/

            return true;

		}

        return false;
    }
    
    void advanceFront (vector<FrontMember>& frontList, vector<Point>& points, double aveTriArea,
                       vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT,
                       PointADT& pointADT, PointADT& edgeCenterADT, EdgeADT& edgeADT, EdgeADT& edge01ADT, int newGridId, vector<Point>& edgeCenters, CircleADT& circleADT, int countr)
    {
        #include "advanceFront.h"
        
        double aveEdgeSize = sqrt ( (4./sqrt(3.) * aveTriArea) );
        //cout << "aveEdgeSize = "  << aveEdgeSize << endl;
        
        // circumradius of average triangle size
        double rho = triEdgeCircumradius (aveEdgeSize, aveEdgeSize, aveEdgeSize);
        //rho *= 1.5;

        int frontCounter = 0;
        while (!frontList.empty())
        {
            bool existingPointPass = false;
            bool circumBoundPass = false;
            bool edgesPass = false;
            bool foundCandPts = false;
            Point CPX;
            int iCPX;
        
            FrontMember& frontFirst = frontList.front();
            int iFrontEdge = frontFirst.edge;
            Edge& frontEdge = edges [ frontFirst.edge ];
            int it0 = frontEdge.t[0];
            int it1 = frontEdge.t[1];
            const Point& t0 = points[ it0 ];
            const Point& t1 = points[ it1 ];

            CVector edgeLen = t1.dim - t0.dim;
            double aveEdgeSize = mag(edgeLen);
            double rho = triEdgeCircumradius (aveEdgeSize, aveEdgeSize, aveEdgeSize);

            cout << "frontList.size() = " << frontList.size() << endl;
            cout << "iFrontEdge = " << iFrontEdge << endl;
            cout << "it0 = " << it0 << endl;
            cout << "it1 = " << it1 << endl;

            //if (countr == 36)
            //{
            //    if (it0 == 34 && it1 == 33)
            //    {
            //        outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
            //        assert(false);
            //    }
            //}

            bool constructedWithExisting = constructWithExisting (foundCandPts, iFrontEdge, existingPointPass, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, edgesPass, t0, t1, circumBoundPass, CPX, iCPX, it0, it1, rho, candPts, frontFirst, frontList, points, edges, triangles, triangleADT, pointADT, edgeADT, edge01ADT, newGridId, circleADT, countr);

            if (!constructedWithExisting)
            //else // if existing point doesn't work
            {
                // get a point normal to front
                Point crP;
                double s = normalDisp (t0.dim, t1.dim, aveTriArea);
                bool gotNewPt = getNewPt (crP, t0, t1, iFrontEdge, s, edge01ADT, triangleADT, points, edges, triangles);
                
                // no need to check circumradius for uniform grid
                // but will be needed in future
                
                // check two forming edges intersections
                //bool newPointPass;
                //if (gotNewPt)
                //{
                //    cout << "checking edges of new point" << endl;
                //    newPointPass = checkTwoFormingEdges (crP, t0, t1, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, edges, edgeADT, points);
                //}

                if (gotNewPt)
                {
                    // add new point to points list
                    addToPointList (crP, points, pointADT);
                    int iCrP = points.size() - 1;

                    int oldiFrontEdge = iFrontEdge;

                    vector <int> edgesAddedToFront;
                    ptCCInter (iFrontEdge, iCrP, newGridId, edgesAddedToFront, circleADT, triangleADT, triangles, frontList, edges, edgeADT, edge01ADT, points, pointADT);

                    //for (int ppp = 0; ppp < points.size(); ++ppp)
                    //{
                    //    auto& pp = points[ppp];
                    //    if ((pp.dim[0] == crP.dim[0] && pp.dim[1] == crP.dim[1]))
                    //    {
                    //        if (ppp != iCrP)
                    //        {
                    //            cout << "pp.id: " << pp.id << endl;
                    //            cout << "crP.id: " << crP.id << endl;
                    //            cout << "it0.id: " << crP.id << endl;
                    //            outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
                    //            assert(false);
                    //        }
                    //    }
                    //}

                    assert(iFrontEdge == oldiFrontEdge);

                    //bool constructedWithExisting = constructWithExisting (foundCandPts, iFrontEdge, existingPointPass, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, edgesPass, t0, t1, circumBoundPass, CPX, iCPX, it0, it1, rho, candPts, frontFirst, frontList, points, edges, triangles, triangleADT, pointADT, edgeADT, edge01ADT, newGridId, circleADT, countr);

                    srchNearbyPts (crP, points, pointADT, candPts, rho);
                    candPts.push_back (iCrP);

                    // construct new triangle with new point and survived edges of dead triangles
                    // survived edges are the ones that newly added to front
                    for (int e: edgesAddedToFront)
                    {
                        int it0 = edges[e].t[0];
                        int it1 = edges[e].t[1];

                        if (Tanemura_Merriam (it0, it1, points, candPts) != iCrP) {continue;}

                        // check edge intersections
                        bool A_CPX_exists, B_CPX_exists;
                        int iA_CPX, iB_CPX;
                        if (checkTwoFormingEdges (crP, points[it0], points[it1], A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, edges, edgeADT, points))
                        {
                            // if pass construct triangle                    
                            construct (iCrP, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, it0, it1, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points, circleADT, e);
                        }
                    }

                    A_CPX_exists = false;
                    B_CPX_exists = false;

                    //if (!constructedWithExisting)
                    {
                        // construct triangle if all previous stages are passed                    
                        construct (iCrP, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, it0,
                                it1, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points, circleADT, iFrontEdge);
                    }
                    //else
                    //{
                    //    erasePoint (iCrP, points, edges, triangles);
                    //    bool success = pointADT.removeViaID (iCrP);
                    //}
                }
                else
                {
                    if (iFrontEdge == 782)
                    {
                        for (auto aa: frontList) {
                            assert(false);
                        }
                    }
                    bool chooseWorstExisPt = false;
                
                    if (foundCandPts)
                    {
                        if (!circumBoundPass)
                        {
                            edgesPass = checkTwoFormingEdges (CPX, t0, t1, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, edges, edgeADT, points);
                            
                            if (edgesPass)
                            {
                                construct (iCPX, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, it0,
                                           it1, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points, circleADT, iFrontEdge);
                                           
                                chooseWorstExisPt = true;
                            }
                        }                        
                    }                    
                
                    if (chooseWorstExisPt == false)
                    {
                        cout << "none of the ways work in AFT::advanceFront(...)" << endl;
                    
                        cout << "it0 = " << it0 << endl;
                        cout << "it1 = " << it1 << endl;
                        cout << "iCPX = " << iCPX << endl;
                        
                        outputTrianglesVTK (points, triangles, "../out", "tri.vtk");

                        exit(-2);
                    }
                cout << "EEEEEEEEEEEEEE" << endl;
                }
            }            
          
            //cout << "frontListSize = " << frontList.size() << endl;
            //cout << "triangles.size() = " << triangles.size() << endl;
            
            //cout << "2rho = "  << 2.*rho << endl;
            /*if (triangles.size() == 140)
            {
                    outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
                    exit(-2);
            }*/
            
            //if (iFrontEdge == 782)
            //{
                        //outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
                        //assert(false);
            //}
            std::string fn = "tri-";
            fn.append(std::to_string(frontCounter));
            fn.append(".vtk");
            outputTrianglesVTK (points, triangles, "../out", fn);
            ++frontCounter;
        }
    }
}
