#ifndef SECULAR_H
#define SECULAR_H
#include "eigenmatrix.h"
class SECU
{
    public:
        nonleaf *Q;
        double  *Lam;
};

SECU *secular(double *d, int dSize, double *v, int vSize,double N);

#endif