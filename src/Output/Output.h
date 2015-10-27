/* 
 * File:   Output.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 6:37 AM
 */

#ifndef OUTPUT_H
#define	OUTPUT_H

#include "../Grid/Grid.h"
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <string>
#include <sys/stat.h>

using std::ofstream;
using std::setw;

string createOutputDir();
template <class T> void log (string directory, T vali, string vals, string unit)
{
    ofstream out;
    string dir = directory;    
    string temps = "log.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);    
    out.open (dir, ofstream::app);
    
    out << vals << " = " << vali;
    
    if (unit.size() != 0)
    {
        out << " "  << unit << endl;
    }
    
    out.close();
}

#endif	/* OUTPUT_H */

