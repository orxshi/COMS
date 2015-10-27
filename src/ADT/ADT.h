/* 
 * File:   ADT.h
 * Author: Orhan Shibliyev
 *
 * Created on July 16, 2014, 6:32 PM
 */

#ifndef ADT_H
#define	ADT_H

#include <vector>
#include <algorithm>
#include <iostream>
#include "../Vector/Vector.h"
#include "../Constants.h"
#include "../Point/Point.h"

#define ADT_DIM 3
#define ADT_VAR (2*ADT_DIM)

using std::vector;
using std::min;
using std::max;
using std::fill;
using std::cout;
using std::endl;
using std::cin;
using std::move;

struct ADT
{
    struct ADTPoint
    {
        int idx;
        vector <CVector> vertices;
        Vector <ADT_VAR> dim;
        /*
         * 0: xmin
         * 1: xmax
         * 2: ymin
         * 3: ymax
         * 4: zmin
         * 5: zmax
         */
    };
    
    int belonging;
    unsigned int sizeTree;
    vector<ADTPoint> points;
    int nIntersections;
    bool searchForNIntersections;
    vector <int> ids;    
    vector <int> idsInTree;
    

    struct Node
    {
        Vector<ADT_VAR> c, d;
        //ADTPoint p;
        ADTPoint* p;
        double key;
        Node* left;
        Node* right;
        unsigned int level;
        bool isEmpty;

        Node()
        {
            left = NULL;
            right = NULL;
            p = NULL;
            isEmpty = true;
        }
    };

    Node *root;
    vector <Node*> searchStack;    
    vector <Node*> addresses;
    vector <Node*> addrsInTree;
    
    //

    void destroy_tree (Node *&leaf);
    void destroy_tree();

    void insert (ADTPoint& point, Node *node, bool& isInserted);

    bool doRegionOverlap (const unsigned int j, const Node *node, const ADTPoint& targetPoint);

    virtual bool compareFunction (const Node *node, const ADTPoint& targetPoint) {}
    virtual bool doCubesOverlap (const Node *node, const ADTPoint& targetPoint) {}

    void searchChildren (const Node *node, const ADTPoint& targetPoint);

    void search (Node* node, const ADTPoint& targetPoint, int& index);
    int search (const ADTPoint& targetPoint);
    
    void removeSearchAddresses ();
    bool removeViaID (int id);

    // Constructor
    ADT();    
    // Move constructor
    ADT (ADT&& other);
    ADT& operator=(ADT&& other);
    // Destructor
    ~ADT();

    void build();
    void moveNodes (Node *&leaf, Node *&otherLeaf);
};

#endif	/* ADT_H */

