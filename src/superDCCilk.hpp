#pragma once
// #include <superDC.h>

// For development (intellisense) purpose only
#if 0
#define OPEN_CILK 1
#include <bits/stdc++.h>
#include <string.h>
#include "BinTree.h"
#include "superDC.h"
#include "divide2.h"
#include "superdcmv_desc.h"
#include "superdcmv_cauchy.h"
#include "eigenmatrix.h"
#include "secular.h"
#include "Generators.h"
#include <sys/time.h>

extern "C"
{
#include <lapacke.h>
#include <lapack.h>
#include <cblas.h>
}

// #include "/opt/intel/oneapi/mkl/2022.0.2/include/mkl_lapacke.h"
using namespace std;
std::pair<double *, nonleaf **> r_RankOneUpdate(double *Lam, int lamSize, std::pair<int, int> zSize, double *Z, nonleaf **n_leaf, int r);
std::pair<double *, double *> computeLeafEig(std::pair<int, int> dSize, double *D, int i);
void Eig_func(int i);

// Global declaration
#include "superDCGlobals.hpp"
#endif

// CILK related includes
#ifdef OPEN_CILK
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
//#include <cilk/reducer_opadd.h>
#else
#define cilk_spawn
#define cilk_sync
#define cilk_for for
#endif

#ifdef OPEN_CILK
void recursiveSuperDCSolver(int node, int nproc, int level, int cutoff){
    vector<int> children = bt->GetChildren(node);
    level++;
    if (children.empty())
    {
        int i = node-1;
        std::pair<double *, double *> E = computeLeafEig(make_pair(resDvd->dSizes[i].first, resDvd->dSizes[i].second), resDvd->D[i], i);

        Lam[i] = E.first;
        LamSizes[i] = resDvd->dSizes[i].first;

        Q0[i] = new EIG_MAT();
        Q0[i]->Q0_leaf = E.second;
        Q0[i]->Q0_nonleaf = NULL;
        q0Sizes[i] = {resDvd->dSizes[i].first, resDvd->dSizes[i].second};

        return;
    }

    int left = children[0];
    int right = children[1];
    
    if (level<cutoff){
    // conquer left child    
    cilk_spawn recursiveSuperDCSolver(left, nproc, level, cutoff);

    // conquer right child
    cilk_spawn recursiveSuperDCSolver(right, nproc, level, cutoff);

    cilk_sync;
    } else{
    // conquer left child    
     recursiveSuperDCSolver(left, nproc, level, cutoff);

    // conquer right child
     recursiveSuperDCSolver(right, nproc, level, cutoff);
    }

    // merge results
    int i = node-1;
    left = children[0];
    right = children[1];
    resDvd->Z[i] = superdcmv_desc(Q0, q0Sizes, (resDvd->Z[i]), resDvd->zSizes[i], bt, i, 1, l, fmmTrigger);
    Lam[i] = new double[(LamSizes[left - 1]) + (LamSizes[right - 1])];
    std::copy(Lam[left - 1], Lam[left - 1] + LamSizes[left - 1], Lam[i]);
    std::copy(Lam[right - 1], Lam[right - 1] + LamSizes[right - 1], Lam[i] + LamSizes[left - 1]);

    LamSizes[i] = (LamSizes[left - 1]) + (LamSizes[right - 1]);
    delete[] Lam[left - 1];
    delete[] Lam[right - 1];

    LamSizes[left - 1] = 0;
    LamSizes[right - 1] = 0;

    int r = resDvd->zSizes[i].second;
    // cout << "R is "<<r<<"\n";

    Q0[i] = new EIG_MAT();
    Q0[i]->Q0_leaf = NULL;

    nonleaf **n_leaf = new nonleaf *[r];
    std::pair<double *, nonleaf **> result = r_RankOneUpdate(Lam[i], LamSizes[i], resDvd->zSizes[i], resDvd->Z[i], n_leaf, r);
    Lam[i] = result.first;
    Q0[i]->Q0_nonleaf = result.second;
    Q0[i]->n_non_leaf = r;
    q0Sizes[i] = {1, r};
    // }

    return;
}

SDC *cilkSuperDC(GEN *A, BinTree *btree, int *m, int mSize, int nProc, int cutoff)
{

    struct timeval timeStart, timeEnd;
    gettimeofday(&timeStart, 0);

    bt = btree;
    cout << "Reached superDC\n";
    // Dividing Stage
    resDvd = divide2(A, bt, m, mSize);

    cout << "Success Divide\n";
    // Conquering stage
    N = bt->GetNumNodes();

    Q0 = new EIG_MAT *[N];
    q0Sizes = new std::pair<int, int>[N];
    Lam = new double *[N];
    LamSizes = new int[N];
    l = new std::pair<int, int>[N];

    for (int k = 0; k < N; k++)
        l[k] = std::make_pair(0, 0);

    l[0] = {0, m[0] - 1};
    int it = 0;
    int lt = 0;

    for (int k = 0; k < N; k++)
    {
        std::vector<int> ch = bt->GetChildren(k + 1);
        if (ch.size() == 0)
        {
            l[k] = {lt, lt + m[it] - 1};
            lt = l[k].second + 1;
            it = it + 1;
        }
        else
        {
            l[k] = {l[ch[0] - 1].first, l[ch[1] - 1].second};
        }
    }

    for (int k = 0; k < N; k++)
        q0Sizes[k] = std::make_pair(0, 0);

    // vector<int> counter(N+1, 0);
    // std::vector<int> WorkQueue(bt->leaves.begin(), bt->leaves.end());

    //  struct timeval timeStart, timeEnd;
    // gettimeofday(&timeStart, 0);
    // cout<<"Number of processors:"<<omp_get_num_procs()<<endl;
    recursiveSuperDCSolver(bt->nodeAtLvl[0][0], nProc, 0, cutoff);

    gettimeofday(&timeEnd, 0);
    long long elapsed = (timeEnd.tv_sec - timeStart.tv_sec) * 1000000LL + timeEnd.tv_usec - timeStart.tv_usec;
    printf("\nDone. %f usecs\n", elapsed / (double)1000000);

    std::ofstream txtOut;
    txtOut.open("output.txt", std::ofstream::out);
    txtOut << setprecision(10) << elapsed / (double)1000000 << " seconds" << endl;
    vector<double> tempeig;
    for (int k = 0; k < LamSizes[N - 1]; k++)
        tempeig.push_back(Lam[N - 1][k]);

    std::sort(tempeig.begin(), tempeig.end());
    int count = 0;
    for (int k = 0; k < LamSizes[N - 1]; k++)
    {
        count++;
        txtOut << setprecision(20) << tempeig[k] << endl;
    }
    // cout << count;
    SDC *resSDC = new SDC();

    resSDC->Q = Q0;
    resSDC->L = Lam[N - 1];
    resSDC->lSize = LamSizes[N - 1];
    resSDC->qSizes = q0Sizes;
    return resSDC;
}
#endif