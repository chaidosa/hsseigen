#include "superdcmv_node.h"
#include "BinTree.h"

void superdcmv_node(double **Q,std::pair<int, int>*qSize, double *X,std::pair<int, int>xSize,BinTree *bt, int index, int ifTrans,std::pair<int,int>* indexRange, double N){
/*
%%% Input:
%%% Qi: hss sturctured cauchylike eigenmatrix of descendants of i
%%% X: vectors to be multiplied
%%% ifTrans = 0: not transpose
%%% ifTrans = 1: transpose
%%% N: size threshold to use fmm
%%% bt: hss tree

%%% Output
%%% ifTrans = 0: X = Qi * X 
%%% ifTrans = 1: X = Qi^T * X;
*/



}