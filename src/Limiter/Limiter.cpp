#include "Limiter.h"

void minMod(const Vector2D<3,N_VAR>& gradL, const Vector2D<3,N_VAR>& gradR, Vector2D<3,N_VAR>& grad)
{
    // http://www.cfdbooks.com/cfdcodes/oned_euler_v1.f90

    double tmp;

    for (int k=0; k<N_VAR; ++k)
    {
        for (int d=0; d<N_DIM; ++d)
        {
            tmp = gradL[k][d] * gradR[k][d];

            if ( tmp > 0. )
            {
                if ( fabs(gradL[k][d]) < fabs(gradR[k][d]) )
                {
                    grad[k][d] = gradL[k][d];
                }
                else
                {
                    grad[k][d] = gradR[k][d];
                }
            }
            else
            {
                grad[k][d] = 0.;
            }
        }
    }
}
