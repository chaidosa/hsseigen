#include "superdcmv_cauchy.h"
#include <iostream>
#include <utility>
#include "eigenmatrix.h"
#include "bsxfun.h"
void superdcmv_cauchy(nonleaf *Q,std::pair<int, int>qSize, double *X,std::pair<int, int>xSize,int ifTrans,double N){
/*
%%% Input:
%%% Q: hss sturctured cauchylike eigenmatrix
%%% X: vectors to be multiplied
%%% ifTrans = 0: not transpose
%%% ifTrans = 1: transpose
%%% N: size threshold to use fmm

%%% Output
%%% ifTrans = 0: X = Q * X 
%%% ifTrans = 1: X = Q^T * X;
*/
    if(ifTrans == 0){

    }
    else{
        // 1st deflation permutation

        //conjugate normalizer
        double *tempX = X + (Q->n1)*xSize.second;
        bsxfun('T',tempX,{xSize.first + Q->n1,xSize.second},Q->v2c,Q->v2cSize);

        //eigenvalue sorting permutation


        //Givens rotation

        //2nd deflation permutation

        //orthogonal Cauchy eigenmatrix
        


    }


}


