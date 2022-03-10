#ifndef SUPERDC_H
#define SUPERDC_H
#include"BinTree.h"
#include "test_mat2hsssym.h"
#include "band2hss.h"
class SDC
{

public:
    double *L, **Q;
    std::pair<int,int> *qSizes;
};

SDC* superDC(tHSSMat *A,  BinTree* bt, int * m, int mSize); 


#endif