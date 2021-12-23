#ifndef SUPERDCMV_DESC_H
#define SUPERDCMV_DESC_H
#include "BinTree.h"
#include "superDC.h"

void superdcmv_desc(double **Q,std::pair<int, int>*qSize, double *X,std::pair<int, int>xSize,BinTree *bt, int index, int ifTrans,std::pair<int,int>* indexRange,double N=1024);

#endif