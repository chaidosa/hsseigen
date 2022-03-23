#ifndef DIVIDE2_H
#define DIVIDE2_H
#include "BinTree.h"
#include "superDC.h"
#include "Generators.h"
class DVD{
    public:
        double **D, **Z;
        std::pair<int,int> *dSizes, *zSizes;

        ~DVD(){
            delete[] D;
            delete[] Z;
            delete[] dSizes;
            delete[] zSizes;
        }

};

DVD* divide2(GEN *A, BinTree *bt,int* m, int mSize);
void norm_svd(double * A, std::pair<int, int>aSize,double norm);
#endif