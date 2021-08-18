#include "Vector.h"

double dotP (const Vector<3>& V1, const Vector<3>& V2)
{
    //////
    double tmp;
    
    tmp = V1[0] * V2[0]
        + V1[1] * V2[1]
        + V1[2] * V2[2];

    return tmp;
}

Vector<3> crossP (const Vector<3>& V1, const Vector<3>& V2)
{
    Vector<3> V3;

    V3[0] = V1[1] * V2[2] - V1[2] * V2[1];
    V3[1] = V1[2] * V2[0] - V1[0] * V2[2];
    V3[2] = V1[0] * V2[1] - V1[1] * V2[0];

    return V3;
}
