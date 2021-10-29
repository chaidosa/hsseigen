#ifndef TEST_MAT2HSSSYM_H
#define TEST_MAT2HSSSYM_H
#include"BinTree.h"
class tHSSMat{
public:
double** D, **U, **R, **B;
std::pair<int,int> *uSizes, *rSizes, *bSizes, *dSizes;
//int N; //Number of elements in {U,R,B}= N-1, {D}=N/2. 
};

tHSSMat* t_mat2hsssym(double* A, int aSize, BinTree* bt, int * m, int mSize, char* tol="tol", double par=0.000001);

#endif