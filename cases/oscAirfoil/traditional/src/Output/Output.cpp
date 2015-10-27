#include "Output.h"

void outLiftCoef (Coeffs& coeffs, double alpha, double time)
{
    if (coeffs.out.is_open())
    {
        coeffs.out << setw(20) << time;
        coeffs.out << setw(20) << alpha;
        coeffs.out << setw(20) << coeffs.liftCoef;
        coeffs.out << endl;
    }
    else
    {
        cout << "coeffs.out is not open in outLiftCoef(...)" << endl;
        exit(-2);
    }
}
