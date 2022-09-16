#ifndef SUPERDCMV_H
#define SUPERDCMV_H
#include<utility>
#include"BinTree.h"
#include"eigenmatrix.h"
double* superdcmv(EIG_MAT **Qt, std::pair<int, int>*qSize, double *x, std::pair<int, int> xSize, BinTree* bt, int * m, int ifTrans,double N);
#endif
