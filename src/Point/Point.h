/* 
 * File:   Point.h
 * Author: orxan
 *
 * Created on August 26, 2014, 8:36 PM
 */

#ifndef POINT_H
#define	POINT_H

#include <forward_list>
#include "../Vector/Vector.h"

using std::forward_list;

struct Point
{
    CVector dim;
    int belonging;
    //bool newlyCreated;
    
    Point();
};

#endif	/* POINT_H */

