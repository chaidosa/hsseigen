#ifndef DIVIDE2_H
#define DIVIDE2_H
#include "test_mat2hsssym.h"
#include "BinTree.h"
#include "superDC.h"

class DVD{
    public:
        double **D, **B, **Z;
        std::pair<int,int> *dSizes, *bSizes, *zSizes;

};

DVD* divide2(tHSSMat *A, BinTree *bt,int* m, int mSize);

#endif