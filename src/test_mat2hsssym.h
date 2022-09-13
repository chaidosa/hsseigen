#ifndef TEST_MAT2HSSSYM_H
#define TEST_MAT2HSSSYM_H
#include"BinTree.h"
#include "Generators.h"
class tHSSMat{
public:
double** D, **U, **R, **B;
std::pair<int,int> *uSizes, *rSizes, *bSizes, *dSizes;
//int N; //Number of elements in {U,R,B}= N-1, {D}=N/2. 
};

GEN* t_mat2hsssym(double* A, int aSize, BinTree* bt, int * m, int mSize, char const* tol="tol", double par=1.0e-10);

#endif
