#include "Init.h"

void OscInit::log (string dir)
{
	ofstream out;
    
    string temps = "log.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);    
    out.open (dir);
    
    out << "rhoInf = " << rhoInf << endl;
	out << "pInf = " << pInf << endl;
	out << "Mach = " << Mach << endl;
	out << "aoa = " << aoa << endl;
	
	out.close();
}

void OscInit::print()
{
	cout << "rhoInf = " << rhoInf << endl;
	cout << "pInf = " << pInf << endl;
	cout << "Mach = " << Mach << endl;
	cout << "aoa = " << aoa << endl;
}
