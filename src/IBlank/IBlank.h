/* 
 * File:   IBlank.h
 * Author: orhan
 *
 * Created on July 18, 2015, 11:27 PM
 */

#ifndef IBLANK_H
#define	IBLANK_H

#include "../Gradient/Gradient.h"

struct Iblank
{
    enum cellCriter_t {WALL=0, SIZE=1};
    cellCriter_t cellCriter;
    
    Iblank ();
    void identify (Grid& grAct, Grid& grPas);
    void interpolate (Grid& gr, Gradient& gradient);
};

#endif	/* IBLANK_H */

