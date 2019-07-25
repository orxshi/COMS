#ifndef GRADIENT_H
#define	GRADIENT_H

#include "../Grid/Grid.h"

struct Gradient
{
    int rank;
    int nProcs;
    int localSize;
    int localSizeNVARNDIM;
    int* localSizes;
    int* localSizesNVARNDIM;
    int* displs;
    int* displsNVARNDIM;
    
    
    vector < array<Vector<N_DIM>, N_VAR> > grad;
    vector <vector<double>> Wx;
    vector <vector<double>> Wy;
    vector <vector<double>> Wz;

    Gradient (Grid& gr);
    
    void initParallelVars (Grid& gr);
    
    void leastSquaresCoeffs (Grid& gr);
    void leastSquaresGrad (Grid& gr);
};

#endif /* GRADIENT_H */
