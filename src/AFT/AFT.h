#ifndef AFT_H
#define	AFT_H

#include <algorithm>
#include <memory>
#include <iomanip>
#include <deque>
#include "../Grid/Grid.h"

using std::sort;
using std::addressof;
using std::cref;
using std::unique_ptr;
using std::setw;
using std::isfinite;
using std::deque;
using std::remove_if;

namespace AFT
{
    struct Point
    {
        bool alive;
        CVector dim;
        int belonging;
        bool newlyCreated;
        vector<int> tri;
        vector<int> e;
        int id;
        
        Point();
        bool eraseParentEdge (int iEdge);
        bool eraseParentTri (int iTri);
    };

    struct Edge
    {
        bool alive;
        vector<int> t;
        int belonging;
        bool newlyCreated;
        vector<int> nei;
        //vector<int> tri;
        int id;

        Edge();
        bool eraseParentTri (int iTri);
    };
    
    struct Triangle
    {   
        bool alive;
        vector<int> p;
        vector<int> e;
        vector <int> nei;
        int id;
        
        Triangle ();
        CVector centroid(const vector<Point>& points);
        double qualityScore (const vector<Point>& points, double aveTriArea, bool verbose, bool& passed);
    };
    
    struct FrontMember
    {
        int edge;
        vector<int> ignore;
        int CPfound;
        bool newPointChecked;
        //int cloPtsMaxSize;
        //deque<int> cloPts;
        
        FrontMember();
    };
    
    struct EdgeADT : public ADT
    {
        virtual bool compareFunction (const Node* node, const ADTPoint& targetPoint);
        virtual bool doCubesOverlap (const Node* node, const ADTPoint& targetPoint);
        void build (const vector<Point>& points, const vector<Edge>& edges);
        ADTPoint createADTPoint (const Point& a, const Point& b);        
    };
    
    struct TriangleADT : public ADT
    {
        virtual bool compareFunction (const Node* node, const ADTPoint& targetPoint);
        virtual bool doCubesOverlap (const Node* node, const ADTPoint& targetPoint);
        void build (const EdgeADT& edgeADT);
        bool localCmpFunc (const ADTPoint& node, const ADTPoint& targetPoint);
        ADTPoint createADTPoint (const Triangle& tri, const vector<Point>& points);
        ADTPoint createADTPoint (const Point& p0, const Point& p1, const Point& p2);
    };
    
    struct PointADT : public ADT
    {
        virtual bool compareFunction (const Node* node, const ADTPoint& targetPoint);
        virtual bool doCubesOverlap (const Node* node, const ADTPoint& targetPoint);        
        void build (const vector<Point>& pointss);
        void build ();
        ADTPoint createADTPoint (const CVector& a, const CVector& b);
    };
    
    struct CircleADT : public ADT
    {
        virtual bool compareFunction (const Node* node, const ADTPoint& targetPoint);
        virtual bool doCubesOverlap (const Node* node, const ADTPoint& targetPoint);        
        void build (const EdgeADT& edgeADT);
        ADTPoint createADTPoint (const Triangle& tri, const vector<Point>& points, int idTri);
    };
    
    // Preset
    extern void setPointsEdges (const vector<Grid>& gr, vector<Point>& points, vector<Edge>& edges, vector<Point>& edgeCenters, int newGridId);
    extern void createFrontList (vector<Edge>& edges, vector<FrontMember>& frontList, vector<Point>& points);
    
    // Front    
    extern void sortFrontList (vector<FrontMember>& frontList, const vector<Point>& points, const vector<Edge>& edges);
    extern void eraseFromFrontList (vector<FrontMember>& frontList, int edge);
    extern void addToFrontList (int edge, vector<FrontMember>& frontList, vector<Point>& points, vector<Edge>& edges);
    extern void eraseExistingEdgeFromFrontList (int ie, vector<FrontMember>& frontList);
    extern void advanceFront (vector<FrontMember>& frontList, vector<Point>& points, double aveTriArea,
                       vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT,
                       PointADT& pointADT, PointADT& edgeCenterADT, EdgeADT& edgeADT, EdgeADT& edge01ADT, int newGridId, vector<Point>& edgeCenters, CircleADT& circleADT);
    
    // Edge
    extern Edge createEdge (int indexA, int indexB, int belonging, bool newlyCreated);
    extern int edgeExists (const int ip1, const int ip2, const vector<Point>& points, const vector<Edge>& edges, bool& exists);
    extern bool checkEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT, const vector<Edge>& edges, const vector<Point>& points, bool& exactMatch, int& result);
    extern int checkNumberOfEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT);    
    extern void addToEdgeList (Edge& edge, int iP1, int iP2, vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points);
    extern void eraseDeadEdges (vector<Edge>& edges, vector<Triangle>& triangles, vector<Point>& points);
    extern bool doIntersect (const CVector& p1, const CVector& q1, const CVector& p2, const CVector& q2, bool& exactMatch);
    extern int orientation (const CVector& p, const CVector& q, const CVector& r);
    extern bool onSegment (const CVector& p, const CVector& q, const CVector& r);
    
    // Triangle
    extern Triangle createTriangle (int e1, int e2, int e3, const vector<Edge>& edges, const vector<Point>& points);
    extern bool triangleIntersect (const Triangle& tri, TriangleADT& triangleADT, const vector<Point>& points);    
    extern bool triangleIntersect (const Point& p0, const Point& p1, const Point& p2, TriangleADT& triangleADT);
    extern bool triangleIntersectMultiple (const Point& p0, const Point& p1, const Point& p2, TriangleADT& triangleADT);
    extern Triangle createTriangle (int e1, int e2, int e3, const vector<Edge>& edges, const vector<Point>& points);
    extern double areaTriangle (const Triangle& tri, const vector<Point>& points);
    extern double charTriangleLength (const Triangle& tri, const vector<Point>& points);    
    extern void triPtsCircums (const CVector& A, const CVector& B, const CVector& C, CVector& cnt, double& radius);
    extern double triEdgeCircumradius (double a, double b, double c);
    extern void flip (vector<Triangle>& triangles, vector<Edge>& edges, const vector<Point>& points);
    extern void outputTrianglesVTK (const vector<Point>& points, const vector<Triangle>& triangles, string dir, string fileName);
    extern void addToTriangleList(vector<Triangle>& triangles, Triangle& tmpTriangle, TriangleADT& triangleADT, vector<Point>& points, CircleADT& circleADT, vector<Edge>& edges, const vector<FrontMember>& frontList);
    extern CVector cntTriangle3Pts (const CVector& p0, const CVector& p1, const CVector& p2);
    extern bool triQuality (const CVector& p0, const CVector& p1, const CVector& p2, double rho);
    extern void eraseDeadTriangles (vector<Triangle>& triangles, vector<Point>& points, vector<Edge>& edges);
    extern double getAveTriArea (const vector<Edge>& edges, const vector<Point>& points);
    extern void outputTriangles (const vector<Point>& points, const vector<Triangle>& triangles);
    
    // Reblanking
    extern void fieldToFringe (Grid& gr, Grid& ogr, int crt);
    extern void fringeToField (Grid& gr, Grid& ogr, int crt);
    
    // Point    
    extern bool getNewPt (Point& crP, const Point& t0, const Point& t1, int iEdge, double aveTriSize, EdgeADT& edge01ADT, TriangleADT& triangleADT, vector<Point>& points, vector<Edge>& edges, vector<Triangle>& triangles);
    extern bool rayCasting (const Point& p, EdgeADT& edgeADT);    
    extern void addToPointList (Point& p, vector<Point>& points, PointADT& pointADT);
    extern double getPointDistance (double r);    
    extern bool pointsNearby (const CVector& range1, const CVector& range2, PointADT& pointADT, PointADT& edgeCenterADT);
    extern bool pointExists (const Point& p, PointADT& pointADT, int& result);
    extern void srchCandPts (FrontMember& fm, vector<Edge>& edges, vector<Point>& points, PointADT& pointADT, deque<int>& candPts, double rho, EdgeADT& edgeADT, EdgeADT& edge01ADT, TriangleADT& triangleADT);
    extern void srchNearbyPts (const Point& p, const vector<Point>& points, PointADT& pointADT, deque<int>& candPts, double rho);
    extern bool checkTwoFormingEdges (const Point& CPX, const Point& A, const Point& B, bool& A_CPX_exists, bool& B_CPX_exists, int& iA_CPX, int& iB_CPX,
            vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points);
    extern bool checkCircumBound (const Point& CPX, const Point& A, const Point& B, double rho);    
    extern void ptCCInter (int baseEdge, const int iCPX, int newGridId, vector <int>& edgesAddedToFront, CircleADT& circleADT, TriangleADT& triangleADT, vector<Triangle>& triangles, vector<FrontMember>& frontList, vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points, PointADT& pointADT);
    extern void eraseDeadPoints (vector<Point>& points, vector<Edge>& edges, vector<Triangle>& triangles);
    
    // New Grid
    extern bool faceExists (const Face& nf, const vector<Face>& face, const vector<::Point>& point, int& index);
    extern void createCells (double offsetZ, const vector<Point>& points, Grid& newGrid, const vector<Triangle>& triangles, int phys, int newGridId);
    
    // Final Grid
    extern bool pointExistsForCreateCells(const ::Point& refPoint, const vector<::Point>& points, int& index);
    extern void createFinalGrid (Grid& finalGrid, const vector<Grid>& gr, const Grid& newGrid);
    
    // AFT
    extern void aft (vector<Grid>& gr, Grid& finalGrid);
    extern void construct (int iCPX, bool A_CPX_exists, bool B_CPX_exists, int iA_CPX, int iB_CPX, int iA, int iB, vector<FrontMember>& frontList,
             vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, vector<Point>& points, CircleADT& circleADT, int iFrontEdge);
    extern void exportToGMSH (const vector<Point>& points, const vector<Edge>& mesh0Edges, const vector<Edge>& mesh1Edges, string dir);
    extern double spacingFnc (double b, double aveTriSize);
    extern int Tanemura_Merriam (int iA, int iB, vector<Point>& points, deque<int> pts);
};

#endif	/* AFT_H */

