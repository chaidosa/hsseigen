#ifndef SUPERDCMV_NODE_H
#define SUPERDCMV_NODE_H
#include "BinTree.h"
#include "eigenmatrix.h"
double *superdcmv_node(EIG_MAT *Qt,std::pair<int, int>qSize, double *tempX,std::pair<int, int>& xSize,BinTree *bt, int index, int ifTrans,double N=17000);

#endif
