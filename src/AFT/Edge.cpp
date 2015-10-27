#include "AFT.h"

namespace AFT
{
    Edge::Edge()
    {
        newlyCreated = false;
        alive = true;
    }
    
    void EdgeADT::build (const vector<Point>& points, const vector<Edge>& edges)
    {
        this->points.resize ( edges.size() );
        
        for (unsigned int e=0; e<edges.size(); ++e)
        {
            const Point& t0 = points[ edges[e].t[0] ];
            const Point& t1 = points[ edges[e].t[1] ];
            
            this->points[e] = this->createADTPoint (t0, t1);
            this->points[e].idx = e;
        }
        
        ADT::build();
    }
    
    ADT::ADTPoint EdgeADT::createADTPoint (const AFT::Point& a, const AFT::Point& b)
    {
        ADTPoint vec;

        for (unsigned int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = min( a.dim[i], b.dim[i] );
            vec.dim[i*2+1] = max( a.dim[i], b.dim[i] );
        }

        vec.vertices.push_back (a.dim);
        vec.vertices.push_back (b.dim);

        return vec;
    }
    
    bool EdgeADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
    {
        bool insideCube = true;

        for (unsigned int d=0; d<ADT_DIM; ++d)
        {
            if (!(node->p->dim[d*2] <= targetPoint.dim[d*2+1]) || !(node->p->dim[d*2+1] >= targetPoint.dim[d*2]))
            {
                insideCube = false;
                break;
            }
        }
    
        return true;
    }
    
    bool EdgeADT::compareFunction (const Node* node, const ADTPoint& targetPoint)
    {
        bool exactMatch;        
        bool inside = doIntersect (targetPoint.vertices[0], targetPoint.vertices[1], node->p->vertices[0], node->p->vertices[1], exactMatch);
        
        return inside;
    }
    
    Edge createEdge (int indexA, int indexB, int belonging, bool newlyCreated)
    {
        Edge edge;
                
        edge.t.push_back (indexA);
        edge.t.push_back (indexB);
        edge.belonging = belonging;
        edge.newlyCreated = newlyCreated;

        return edge;
    }
    
    int edgeExists (const int ip1, const int ip2, const vector<Point>& points, const vector<Edge>& edges, bool& exists)
    {
        exists = false;
        
        for (unsigned int e=0; e<edges.size(); ++e)
        {
            int it0 = edges[e].t[0];
            int it1 = edges[e].t[1];
            
            if ( (it0 == ip1) || (it1 == ip1) )
            {
                if ( (it0 == ip2) || (it1 == ip2) )
                {
                    exists = true;
                    return e;
                }
            }
        }
    }
    
    bool checkEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT, const vector<Edge>& edges, const vector<Point>& points, bool& exactMatch, int& result)
    {
        result = -1;

        ADT::ADTPoint vec = edgeADT.createADTPoint (frontListPoint, closestPoint);

        edgeADT.searchForNIntersections = false;
        result = edgeADT.search (vec);

        exactMatch = false;
        if (result != -1)
        {
            //cout << result << endl;
            //cout << edges.size() << endl;
            doIntersect (frontListPoint.dim, closestPoint.dim, points[edges[result].t[0]].dim, points[edges[result].t[1]].dim, exactMatch);
            //cout << "result = " << result << endl;
            //cout << "edges[result].alive = " << edges[result].alive << endl;
            //cout << "edges[result].t[0] = " << edges[result].t[0] << endl;
            //cout << "edges[result].t[1] = " << edges[result].t[1] << endl;
            return true;
        }
        
        return false;
    }
    
    int checkNumberOfEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT)
    {
        ADT::ADTPoint vec = edgeADT.createADTPoint (frontListPoint, closestPoint);

        edgeADT.searchForNIntersections = true;
        edgeADT.search (vec);
        
        return edgeADT.nIntersections;
    }
    
    void addToEdgeList (Edge& edge, int iP1, int iP2, vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points)
    {
        bool tempBool;
        
        edges.push_back (edge);
        ADT::ADTPoint vec = edgeADT.createADTPoint (points[iP1], points[iP2]);
        vec.idx = edges.size() - 1;
        
        edgeADT.insert (vec, edgeADT.root, tempBool);
        
        points[iP1].e.push_back (edges.size() - 1);
        points[iP2].e.push_back (edges.size() - 1);
    }
    
    void eraseDeadEdges (vector<Edge>& edges, vector<Triangle>& triangles, vector<Point>& points)
    {
        for (int ie=0; ie<edges.size(); ++ie)
        {
            Edge& e = edges[ie];
            
            e.id = ie;
        }
        
        edges.erase (remove_if (edges.begin(), edges.end(), [](Edge& e) { return e.alive == false; }), edges.end());
        
        /*for (int ie=0; ie<edges.size(); ++ie)
        {
            Edge& e = edges[ie];
        
            if (e.alive == false)
            {
                edges.erase(edges.begin() + ie);
            }
        }*/
        
        for (int ie=0; ie<edges.size(); ++ie)
        {
            Edge& e = edges[ie];
            
            for (int ip: e.t)
            {
                Point& p = points[ip];
                
                for (int pp=0; pp<p.e.size(); ++pp)
                {
                    if (p.e[pp] == e.id)
                    {
                        p.e[pp] = ie;
                        break;
                    }
                }
            }
            
            for (int t: e.nei)
            {            
                if (t != -1)
                {
                    Triangle& tri = triangles[t];
                    
                    for (int pp=0; pp<tri.p.size(); ++pp)
                    {
                        if (tri.e[pp] == e.id)
                        {
                            tri.e[pp] = ie;
                            break;
                        }
                    }
                }
            }
        }
        
        for (int ie=0; ie<edges.size(); ++ie)
        {
            Edge& e = edges[ie];
            
            e.id = ie;
        }
    }
    
    bool Edge::eraseParentTri (int iTri)
    {
        for (int it=0; it<nei.size(); ++it)        
        {
            int t = nei[it];
        
            if (t == iTri)
            {
                nei.erase (nei.begin() + it);
                return true;
            }
        }
        
        return false;
    }
    
    bool doIntersect (const CVector& p1, const CVector& q1, const CVector& p2, const CVector& q2, bool& exactMatch)
    {
        // http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/

        // The main function that returns true if line segment 'p1q1'
        // and 'p2q2' intersect

        // Find the four orientations needed for general and
        // special cases

        int o1 = orientation(p1, q1, p2);
        int o2 = orientation(p1, q1, q2);
        int o3 = orientation(p2, q2, p1);
        int o4 = orientation(p2, q2, q1);    
        
        

        /*if (o1 == -1 || o2 == -1 || o3 == -1 || o4 == -1)
        {
            return true;
        }*/
        // General case
        if (o1 == 0 || o2 == 0 || o3 == 0 || o4 == 0)
        {
            // Special Cases
            // p1, q1 and p2 are colinear and p2 lies on segment p1q1
            if (o1 == 0 && onSegment(p1, p2, q1)) return true;

            // p1, q1 and p2 are colinear and q2 lies on segment p1q1
            if (o2 == 0 && onSegment(p1, q2, q1)) return true;

            // p2, q2 and p1 are colinear and p1 lies on segment p2q2
            if (o3 == 0 && onSegment(p2, p1, q2)) return true;

             // p2, q2 and q1 are colinear and q1 lies on segment p2q2
            if (o4 == 0 && onSegment(p2, q1, q2)) return true;

            unsigned int count = 0;
            exactMatch = false;

            if ( onSegment(p1, p2, q1) ) ++count;
            if ( onSegment(p1, q2, q1) ) ++count;
            if ( onSegment(p2, p1, q2) ) ++count;
            if ( onSegment(p2, q1, q2) ) ++count;

            if (count == 0)
            {
                if ( ((p1[0] == p2[0]) && (p1[1] == p2[1])) || ((p1[0] == q2[0]) && (p1[1] == q2[1])) )
                {
                    if ( ((q1[0] == p2[0]) && (q1[1] == p2[1])) || ((q1[0] == q2[0]) && (q1[1] == q2[1])) )
                    {
                        exactMatch = true;
                        return true;
                    }
                }
            }
        }
        else if (o1 != o2 && o3 != o4)
        {
            //if (o1*o2*o3*o4 != 0)
            //{
                return true;
            //}
        }

        return false; // Doesn't fall in any of the above cases
    }
    
    int orientation (const CVector& p, const CVector& q, const CVector& r)
    {
        // To find orientation of ordered triplet (p, q, r).
        // The function returns following values
        // 0 --> p, q and r are colinear
        // 1 --> Clockwise
        // 2 --> Counterclockwise

        // See 10th slides from following link for derivation of the formula
        // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf

        double val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);

        if (val == 0.) return 0;  // colinear

        return (val > 0.) ? 1 : 2; // clock or counterclock wise
    }
    
    bool onSegment (const CVector& p, const CVector& q, const CVector& r)
    {
        // Given three colinear Vectors p, q, r, the function checks if
        // Vector q lies on line segment 'pr'

        bool AllEqualX;
        bool AllEqualY;
        
        if ( q[0] == p[0] && q[1] == p[1] ) {return false;}
        if ( q[0] == r[0] && q[1] == r[1] ) {return false;}

        if ( q[0] == max(p[0], r[0]) && q[0] == min(p[0], r[0]) )
        {
            AllEqualX = true;
        }
        else
        {
            AllEqualX = false;
        }

        if ( q[1] == max(p[1], r[1]) && q[1] == min(p[1], r[1]) )
        {
            AllEqualY = true;
        }
        else
        {
            AllEqualY = false;
        }

        if ( (((q[0] < max(p[0], r[0])) && (q[0] > min(p[0], r[0]))) || AllEqualX) &&
             (((q[1] < max(p[1], r[1])) && (q[1] > min(p[1], r[1]))) || AllEqualY) )
           return true;

        return false;
    }
}
