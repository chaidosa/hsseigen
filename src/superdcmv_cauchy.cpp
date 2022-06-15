#include "superdcmv_cauchy.h"
#include <iostream>
#include <utility>
#include <string.h>
#include "eigenmatrix.h"
#include "bsxfun.h"
#include "cauchylikematvec.h"

extern "C"{

    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>

}


double * superdcmv_cauchy(nonleaf *Qq,std::pair<int, int>qSize, double *Xx,std::pair<int, int>xSize,int ifTrans,double N){
/*
%%% Input:
%%% Q: hss sturctured cauchylike eigenmatrix
%%% X: vectors to be multiplied
%%% ifTrans = 0: not transpose
%%% ifTrans = 1: transpose
%%% N: size threshold to use fmm

%%% Output:
%%% ifTrans = 0: X = Q * X 
%%% ifTrans = 1: X = Q^T * X;
*/
    nonleaf *Q = Qq;
    double *X = Xx;
    if(ifTrans == 0){
        assert(false);
    }
    
    else{
        // 1st deflation permutation
        int *temp = new int[Q->TSize.first];
        for(int i = 0; i < Q->TSize.first; i++)
            temp[i] = (Q->T[i]) + 1;        
       
        arrange_elements2(X, xSize, (temp), Q->TSize);
        delete [] temp;      

        //double *tempX = X + (Q->n1)*xSize.second;
        double *tempXX = new double[(xSize.first - Q->n1) * xSize.second];
        memcpy(tempXX, X+(Q->n1 * xSize.second), sizeof(double)*((xSize.first - Q->n1) * xSize.second));

        bsxfun('T',&tempXX,{xSize.first - Q->n1,xSize.second},Q->v2c,Q->v2cSize);
       
        //eigenvalue sorting permutation
        int *tempI = new int[Q->ISize.first];        
        for(int row = 0; row < Q->ISize.first; row++)
            tempI[row] = Q->I[row] + 1;

        arrange_elements2(tempXX, {xSize.first - Q->n1, xSize.second}, (tempI), Q->ISize);

        delete[] tempI;
        // Givens rotation
        
        double *tempArr = new double[2*xSize.second];
        double *MulR = new double[2*xSize.second];
        for(int l = 0; l < (Q->GSize.first); l = l + 4){
            int p =(int)Q->G[l];
            int j =(int)Q->G[l+1];
            double c = Q->G[l+2];
            double s = Q->G[l+3];            
           
           for(int cpy = 0; cpy < 2; cpy++){
               memcpy(tempArr+cpy*xSize.second, tempXX + (p + cpy)*xSize.second, sizeof(double)*(xSize.second));
           }
           double GivensArr[2*2] = {c, -s, s, c};
           
           cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, xSize.second, 2, 1, GivensArr, 2, tempArr, xSize.second, 0.0, MulR, xSize.second);
            
           for(int cpy = 0; cpy < 2; cpy++){
               memcpy(tempXX + (p + cpy)*xSize.second, MulR + cpy*(xSize.second), sizeof(double)*(xSize.second));
           }
          
        }

        delete [] tempArr;
        delete [] MulR;
 
        // second deflation permutation
        int *tempJ = new int[Q->JSize.first];
        for(int i = 0; i < Q->JSize.first; i++)
            tempJ[i] = Q->J[i] + 1;

        arrange_elements2(tempXX, {xSize.first - Q->n1, xSize.second}, (tempJ), Q->JSize);
        
        memcpy(X+(Q->n1 * xSize.second), tempXX, sizeof(double)*((xSize.first - Q->n1) * xSize.second));
        delete[] tempXX;
        
        int rowSize = (Q->n)-(Q->n1 + Q->n2 + 1) + 1;
        double **Qc = Q->QC;        
        tempXX = new double[rowSize * xSize.second];
        
        memcpy(tempXX, X+(Q->n1 + Q->n2)*xSize.second, sizeof(double)*(rowSize * xSize.second));
        tempXX = cauchylikematvec(Qc, (Q->qcSizes), Q->Org, (tempXX), {rowSize, xSize.second}, 1, N);
        memcpy(X+(Q->n1 + Q->n2)*xSize.second, tempXX, sizeof(double)*(rowSize * xSize.second));
        delete[] tempXX;
       
    }
    
return X;
}


