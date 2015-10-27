/* 
 * File:   Init.h
 * Author: orhan
 *
 * Created on June 18, 2015, 9:09 PM
 */

#ifndef INIT_H
#define	INIT_H

#include "../Grid/Grid.h"

struct Init
{
    virtual void read() = 0;
    virtual void init (Grid& gr) = 0;
};

#endif	/* INIT_H */

