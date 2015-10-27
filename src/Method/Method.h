/* 
 * File:   Method.h
 * Author: Orhan Shibliyev
 *
 * Created on July 15, 2014, 9:52 PM
 */

#ifndef METHOD_H
#define	METHOD_H

#include <iostream>
#include <fstream>
#include <vector>
#include "../Output/Output.h"
#include "../Time/Time.h"
#include "../AFT/AFT.h"
#include "../Solver/Solver.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::vector;

enum class method_t {UNDEFINED=-1, SINGLE=0, MULTI=1, NEW=2};

void SSG (Grid& gr);

#endif	/* METHOD_H */

