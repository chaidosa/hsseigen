#include "superdcmv_cauchy.h"
#include <iostream>
#include <utility>
#include "eigenmatrix.h"
void superdcmv_cauchy(nonleaf *Q,std::pair<int, int>*qSize, double *X,std::pair<int, int>xSize,int ifTrans,double N){
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



}


