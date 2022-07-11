#ifndef SUPERDC_H
#define SUPERDC_H
#include"BinTree.h"
#include "test_mat2hsssym.h"
#include "band2hss.h"
#include "eigenmatrix.h"
class SDC
{

public:
    double *L;
    EIG_MAT ** Q;
    std::pair<int,int> *qSizes;

    ~SDC(){
        delete [] L;
        delete [] Q;
        delete [] qSizes;
        cout<<"\ndeleted SDC\n"; 
    }
};

extern int fmmTrigger;

SDC* superDC(GEN *A, BinTree* bt, int * m, int mSize, int nProc); 

#endif
