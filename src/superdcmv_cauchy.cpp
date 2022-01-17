#include "superdcmv_cauchy.h"
#include <iostream>
#include <utility>
#include <string.h>
#include "eigenmatrix.h"
#include "bsxfun.h"
#include "cauchylikematvec.h"
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
/*    
    else{
        // 1st deflation permutation
        double *temp;
        arrange_elements2(X,xSize,Q->T,Q->TSize,temp);
        delete [] X;
        X = temp;
        temp = NULL;
        xSize = {Q->TSize.second,xSize.second};
        //conjugate normalizer
        double *tempX = X + (Q->n1)*xSize.second;
        bsxfun('T',&tempX,{xSize.first + Q->n1,xSize.second},Q->v2c,Q->v2cSize);

        //eigenvalue sorting permutation
        double *tempI = Q->I+(Q->n1)*Q->ISize.second;
        double *result_eigen_permutation;
        arrange_elements2(tempX,{xSize.first + Q->n1,xSize.second},tempI,{(Q->n1),Q->ISize.second},result_eigen_permutation);
        memcpy(tempX,result_eigen_permutation,sizeof(double)*((xSize.first + Q->n1)*xSize.second));
        //Givens rotation

        //2nd deflation permutation
        double *res_2_deflation;
        double *tempJ = Q->J+(Q->n1)*Q->JSize.second;
        arrange_elements2(tempX,{xSize.first + Q->n1,xSize.second},tempJ,{Q->n1, Q->JSize.second},res_2_deflation);
        memcpy(tempX,res_2_deflation,sizeof(double)*((xSize.first + Q->n1)*xSize.second));
        //orthogonal Cauchy eigenmatrix
        


    }
    Todo:: modify non-working method
*/

}


