#include "Init.h"

void OscInit::read()
{
	string tmps;
	ifstream in;
    in.open("Init/init.dat");
    
    if (!in.is_open())
    {
    	cout << "Could not open init file in InitData::read()" << endl;
    	exit(-2);
    }
    
    in >> tmps; in >> rhoInf;
    in >> tmps; in >> pInf;
    in >> tmps; in >> Mach;
    in >> tmps; in >> aoa;
    
    in.close();
}

    
