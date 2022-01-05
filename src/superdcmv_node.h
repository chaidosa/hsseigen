#ifndef SUPERDCMV_NODE_H
#define SUPERDCMV_NODE_H
#include "BinTree.h"
#include "eigenmatrix.h"
void superdcmv_node(EIG_MAT *Q,std::pair<int, int>qSize, double **tempX,std::pair<int, int>xSize,BinTree *bt, int index, int ifTrans,double N=1024);

#endif