#ifndef SUPERDC_H
#define SUPERDC_H
#include"BinTree.h"
#include "test_mat2hsssym.h"
#include "band2hss.h"
#include "eigenmatrix.h"

// #define DIST 1 // For development
// #define HYBRD 1 // For Development only
// #define PARALLEL_TASKS 1

class SDC
{

public:
    double *L;
    int lSize;
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

#ifdef DIST
#include <mpi.h>
SDC *dsuperDc(GEN *A, BinTree *btree, int *m, int mSize, int nProc, MPI_Comm process_grid);
#endif

#ifdef HYBRD
#include <mpi.h>
SDC* HybridSuperDC(GEN *A, BinTree* bt, int * m, int mSize, int nProc, MPI_Comm process_grid);
#endif

#ifdef PARALLEL_TASK
SDC *taskSuperDC(GEN *A, BinTree *btree, int *m, int mSize, int nProc);
#endif

#endif
