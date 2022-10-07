#ifndef SUPERDCMV_H
#define SUPERDCMV_H
#include<utility>
#include"BinTree.h"
#include"eigenmatrix.h"
double* superdcmv(EIG_MAT **Qt, std::pair<int, int>*qSize, double *x, std::pair<int, int>* xSize, int* m, int mSize, BinTree* bt, int k, int ifTrans,double N);
#endif
