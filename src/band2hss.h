#ifndef BAND2HSS_H
#define BAND2HSS_H
#include "BinTree.h"
#include "Generators.h"
class B2HSS{
public:
double** D, **U, **R, **B;
std::pair<int,int> *uSizes, *rSizes, *bSizes, *dSizes;
//int N; //Number of elements in {U,R,B}= N-1, {D}=N/2. 
};

GEN *band2hss(double *A, int aSize, BinTree* bt, int* m, int mSize, int w);

#endif