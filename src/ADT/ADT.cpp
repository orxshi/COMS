#include "ADT.h"

ADT::ADT()
{
    sizeTree = 0;
    root = NULL;
    nIntersections = 0;
    searchForNIntersections = false;
}

/*ADT::ADT (ADT&& other)
{
    sizeTree = other.sizeTree;
    //points.clear();
    points = move(other.points);
    //idsInTree.clear();
    idsInTree = move(other.idsInTree);
    
    for (unsigned int i=0; i<other.searchStack.size(); ++i)
    {
        other.searchStack[i] = NULL;
    }
    other.searchStack.clear();
    
    for (unsigned int i=0; i<other.addresses.size(); ++i)
    {
        other.addresses[i] = NULL;
    }
    other.addresses.clear();
    
    //addrsInTree.clear();
    addrsInTree.reserve (other.addrsInTree.size());
    for (unsigned int i=0; i<other.addrsInTree.size(); ++i)
    {        
        addrsInTree.push_back (other.addrsInTree[i]);
        other.addrsInTree[i] = NULL;
    }
    other.addrsInTree.clear();
    
    //moveNodes (root, other.root);
    root = other.root;
    other.root = NULL;
}*/

/*ADT& ADT::operator=(ADT&& other)
{
    sizeTree = other.sizeTree;
    //points.clear();
    points = move(other.points);
    //idsInTree.clear();
    idsInTree = move(other.idsInTree);
    
    for (unsigned int i=0; i<other.searchStack.size(); ++i)
    {
        other.searchStack[i] = NULL;
    }
    other.searchStack.clear();
    
    for (unsigned int i=0; i<other.addresses.size(); ++i)
    {
        other.addresses[i] = NULL;
    }
    other.addresses.clear();
    
    //addrsInTree.clear();
    addrsInTree.reserve (other.addrsInTree.size());
    for (unsigned int i=0; i<other.addrsInTree.size(); ++i)
    {        
        addrsInTree.push_back (other.addrsInTree[i]);
        other.addrsInTree[i] = NULL;
    }
    other.addrsInTree.clear();
    
    //moveNodes (root, other.root);
    root = other.root;
    other.root = NULL;
}*/

void ADT::moveNodes (Node *&leaf, Node *&otherLeaf)
{
    if (otherLeaf != NULL)
    {
        leaf = otherLeaf;
        
        if (otherLeaf->left != NULL)
        {
            moveNodes (leaf->left, otherLeaf->left);
        }
        if (otherLeaf->right != NULL)
        {
            moveNodes (leaf->right, otherLeaf->right);
        }
        
        //leaf = move(otherLeaf);
        otherLeaf = NULL;
    }
}

ADT::~ADT()
{    
    destroy_tree();
}

void ADT::destroy_tree(Node *&leaf)
{
    if (leaf != NULL)
    {
        destroy_tree(leaf->left);
        destroy_tree(leaf->right);
        
        if (leaf->p != NULL)
        {
            delete leaf->p;
        }
        
        delete leaf;
        leaf = NULL;
    }
}

void ADT::destroy_tree()
{   
    for (unsigned int i=0; i<searchStack.size(); ++i)
    {
        searchStack[i] = NULL;
    }
    
    for (unsigned int i=0; i<addrsInTree.size(); ++i)
    {
        addrsInTree[i] = NULL;
    }
    
    for (unsigned int i=0; i<addresses.size(); ++i)
    {
        addresses[i] = NULL;
    }

    destroy_tree(root);
}

void ADT::insert (ADTPoint& point, Node* node, bool& isInserted)
{
    isInserted = false;

    unsigned int j = node->level % ADT_VAR;
    node->key = 0.5 * (node->c[j] + node->d[j]);

    // point is in left region
    if (point.dim[j] < node->key)
    {
        if (node->left == NULL)
        {
            node->left = new Node();
            ++sizeTree;

            for (int d=0; d<ADT_VAR; ++d)
            {
                if (d != j)
                {
                    node->left->c[d] = node->c[d];
                    node->left->d[d] = node->d[d];
                }
                else
                {
                    node->left->c[d] = node->c[d];
                    node->left->d[d] = node->key;
                }
            }

            node->left->level = node->level + 1;
            isInserted = true;
            node->left->p = new ADTPoint (point);
            idsInTree.push_back (point.idx);
            addrsInTree.push_back (node->left);            
            node->left->isEmpty = false;
        }
        else if (node->left->isEmpty)
        {
            ++sizeTree;
            node->left->p = new ADTPoint (point);
            node->left->isEmpty = false;
            idsInTree.push_back (point.idx);
            addrsInTree.push_back (node->left);
            isInserted = true;
        }
        else
        {
            insert(point, node->left, isInserted);
        }
    }
    // point is in right region
    else if (point.dim[j] >= node->key)
    {
        if (node->right == NULL)
        {
            node->right = new Node();
            ++sizeTree;

            for (int d=0; d<ADT_VAR; ++d)
            {
                if (d != j)
                {
                    node->right->c[d] = node->c[d];
                    node->right->d[d] = node->d[d];
                }
                else
                {
                    node->right->c[d] = node->key;
                    node->right->d[d] = node->d[d];
                }
            }

            node->right->level = node->level + 1;
            isInserted = true;
            node->right->p = new ADTPoint (point);
            idsInTree.push_back (point.idx);
            addrsInTree.push_back (node->right);
            node->right->isEmpty = false;
        }
        else if (node->right->isEmpty)
        {
            ++sizeTree;
            node->right->p = new ADTPoint (point);
            node->right->isEmpty = false;
            idsInTree.push_back (point.idx);
            addrsInTree.push_back (node->right);
            isInserted = true;
        }
        else
        {
            insert(point, node->right, isInserted);
        }
    }
    else
    {
        cout << "point.dim[j] = " << point.dim[j] << endl;
        cout << "node->key = " << node->key << endl;
        exit(-2);
    }
}



bool ADT::doRegionOverlap (const unsigned int j, const Node* node, const ADTPoint& targetPoint)
{
    bool regionOverlap = true;

    if ( (j % 2 == 0) )
    {
        if ( !(node->c[j] <= targetPoint.dim[j+1]) )
        {
            regionOverlap = false;
        }
    }
    else
    {
        if ( !(node->d[j] >= targetPoint.dim[j-1]) )
        {
            regionOverlap = false;
        }
    }

    return regionOverlap;
}

void ADT::searchChildren (const Node* node, const ADTPoint& targetPoint)
{
    unsigned int j = node->level % ADT_VAR;

    // check left
    if (node->left != NULL)
    {
        if ( doRegionOverlap (j, node->left, targetPoint) )
        {
            searchStack.push_back (node->left);
        }
    }

    // check right
    if (node->right != NULL)
    {
        if ( doRegionOverlap (j, node->right, targetPoint) )
        {
            searchStack.push_back (node->right);
        }
    }
}

void ADT::search (Node* node, const ADTPoint& targetPoint, int& index)
{
    

    if (node != NULL)
    {
        /*if (targetPoint.idx == 231)
        {
            cout << "alibaba" << endl;
            //exit(-2);
        }*/
    
        // check whether the point is inside the element
        if (!node->isEmpty && node->p->idx!=-1 && doCubesOverlap (node, targetPoint) && compareFunction (node, targetPoint) )
        {
            if (searchForNIntersections)
            {
                ++nIntersections;
                ids.push_back (node->p->idx);
                addresses.push_back (node);
                
                searchChildren (node, targetPoint);

                if (!searchStack.empty())
                {
                    Node* last = searchStack.back();                
                    searchStack.back() = NULL;                
                    searchStack.pop_back();
                    search (last, targetPoint, index);
                }
            }
            else
            {
                fill (searchStack.begin(), searchStack.end(), nullptr);
                searchStack.clear();
            
                index = node->p->idx;
                addresses.push_back (node);
            
                if (node != root)
                {
                    node = NULL;
                }
            }
        }
        else
        {
            searchChildren (node, targetPoint);

            if (!searchStack.empty())
            {
                Node* last = searchStack.back();                
                searchStack.back() = NULL;                
                searchStack.pop_back();
                search (last, targetPoint, index);
            }
        }
    }
}

int ADT::search (const ADTPoint& targetPoint)
{
    int i = -1;
    nIntersections = 0;
    ids.clear();    
    for (int j=0; j<addresses.size(); ++j) { addresses[j] = NULL; } addresses.clear();
    
    fill (searchStack.begin(), searchStack.end(), nullptr);
    searchStack.clear();

    bool regionOverlap = true;

    for (int d=0; d<ADT_DIM; ++d)
    {
        if ( !(root->c[2*d] <= targetPoint.dim[2*d+1]) || !(root->d[2*d+1] >= targetPoint.dim[2*d]) )
        {
            regionOverlap = false;
        }
    }

    if (regionOverlap)
    {
        
        search (root, targetPoint, i);
    }    

    fill (searchStack.begin(), searchStack.end(), nullptr);
    searchStack.clear();    
    
    return i;
}

void ADT::removeSearchAddresses ()
{
    for (int j=0; j<addresses.size(); ++j)
    {
        if (addresses[j] != NULL)
        {
            addresses[j]->isEmpty = true;
            delete addresses[j]->p;
            addresses[j] = NULL;
        }
    }
}

bool ADT::removeViaID (int id)
{
    for (int i=0; i<idsInTree.size(); ++i)
    {
        if (id == idsInTree[i])
        {
            if (addrsInTree[i] != NULL)
            {
                addrsInTree[i]->isEmpty = true;
                
                if (addrsInTree[i]->p != NULL)
                {
                    delete addrsInTree[i]->p;
                    addrsInTree[i]->p = NULL;
                }
                
                addrsInTree[i] = NULL;
                addrsInTree.erase (addrsInTree.begin() + i);
                idsInTree.erase (idsInTree.begin() + i);
                return true;
            }
        }
    }
    
    return false;
}

void ADT::build ()
{
    bool isInserted;
    
    if (root == NULL)
    {
        root = new Node();
        ++sizeTree;
        root->level = 0;

        for (int d=0; d<ADT_VAR; ++d)
        {
            root->c[d] = BIG_POS_NUM;
            root->d[d] = BIG_NEG_NUM;
        }

        for (unsigned int p=0; p<points.size(); ++p)
        {
            for (int d=0; d<ADT_VAR; ++d)
            {
                root->c[d] = min (points[p].dim[d], root->c[d]);
                root->d[d] = max (points[p].dim[d], root->d[d]);
            }
        }

        root->p = new ADTPoint (points.front());        
        root->isEmpty = false;
        idsInTree.push_back (points.front().idx);
        addrsInTree.push_back (root);

        points.erase (points.begin());
    }

    while (!points.empty())
    {
        isInserted = false;
        insert (points.front(), root, isInserted);
        if (isInserted)
        {
            points.erase (points.begin());
        }
        else
        {
            cout << "not inserted" << endl;
            exit(-2);
        }
    }
}
