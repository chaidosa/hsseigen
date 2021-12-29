#include "cauchylikematvec.h"

void cauchylikematvec(double **Qc, std::pair<int,int>*qcSizes,double *X, std::pair<int,int>xSize,int ifTrans, double N=1024){
/*
    %%% Input:
    %%% Q: {v,s,d,lam} stores the cauchylike matrix  v_i * s_j / (d_i - lam_j) 
    %%% X: vectors to be multiplied
    %%% ifTrans = 0: not transpose
    %%% ifTrans = 1: transpose
    %%% N: size threshold to use fmm

    %%% Output
    %%% t = 0: Y = Q * X 
    %%% t = 1: Y = Q^T * X;



*/  
}