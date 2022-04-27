#include "superdcmv_node.h"
#include "BinTree.h"
#include<string.h>
#include <assert.h>
#include "eigenmatrix.h"
#include "superdcmv_cauchy.h"
#ifndef OPENBLAS 
extern "C"
{
#endif
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
#ifndef OPENBLAS 
}
#endif

double *superdcmv_node(EIG_MAT *Qt,std::pair<int, int>qSize, double *tempX,std::pair<int, int>xSize,BinTree *bt, int index, int ifTrans,double N){
/*
%%% Input:
%%% Q: hss sturctured cauchylike eigenmatrix of descendants of i
%%% X: vectors to be multiplied
%%% ifTrans = 0: not transpose
%%% ifTrans = 1: transpose
%%% N: size threshold to use fmm
%%% bt: hss tree

%%% Output
%%% ifTrans = 0: X = Qi * X 
%%% ifTrans = 1: X = Qi^T * X;
*/
    EIG_MAT *Q = Qt;
    double *X = tempX;
    double alpha = 1;
    double beta  = 0;
    std::vector<int> ch = bt->GetChildren(index+1);
    if(ch.size()==0)
    {
        if(ifTrans == 0){
            if(qSize.second != xSize.first){
                cout<<"Multiplication is not possible";
                assert(false);
            }

            double *tempQX = new double[qSize.first*xSize.second];
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,qSize.first,xSize.second,qSize.second,alpha,Q->Q0_leaf,qSize.second,X,xSize.second,beta,tempQX,xSize.second); //Nikhil: memory for Q0_leaf seems to have not been allocated (multiple places). Confirm.
            //delete [] X;
            tempX = tempQX;
            tempQX = NULL;
        }

        else
        {
            if(qSize.first != xSize.first){
                cout<<"Multiplication is not possible";
                assert(false);
            }

            double *tempQX = new double[qSize.second*xSize.second];
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,qSize.second,xSize.second,qSize.first,alpha,Q->Q0_leaf,qSize.second,X,xSize.second,beta,tempQX,xSize.second);
            //delete [] X;
            tempX = tempQX;
            tempQX = NULL;
        }
    }

    else
    {
        int r = qSize.second;
        if(ifTrans == 0)
        {
            //for j = r:-1:1
            //    [X, nflops1] = superdcmv_cauchy(Qi{j}, X, t, N);
        }
        else
        {
            for(int j = 0; j < r; j++)
            {
              double *tempXX = new double[xSize.first * xSize.second];
              memcpy(tempXX, X, sizeof(double)*(xSize.first * xSize.second));
              superdcmv_cauchy(Q->Q0_nonleaf[j],{7 , 1}, tempXX, xSize, 1, 1024);
              memcpy(X, tempXX, sizeof(double)*(xSize.first * xSize.second));  
            }
            tempX = X;
        }
    }
    return tempX;
}
