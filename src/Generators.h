#ifndef GENERATORS_H
#define GENERATORS_H
#include"BinTree.h"
class GEN{
public:
double** D, **U, **R, **B;
std::pair<int,int> *uSizes, *rSizes, *bSizes, *dSizes;
//int N; //Number of elements in {U,R,B}= N-1, {D}=N/2. 
};

GEN* HssGenerators(double *A, int aSize, BinTree* bt, int* m, int mSize, int w = 0, int MorB = 1);
GEN* InitGen(int N);

// #define DIST 1
#if defined(DIST) || defined(HYBRD)
#include <mpi.h>
void sendGen(GEN* send, MPI_Comm process_grid, int N);

#endif

#endif