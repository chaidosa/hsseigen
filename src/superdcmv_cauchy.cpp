#include "superdcmv_cauchy.h"
#include <iostream>
#include <utility>
#include <string.h>
#include "eigenmatrix.h"
#include "bsxfun.h"
#include "cauchylikematvec.h"

/*extern "C"{

    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>

}
*/

void superdcmv_cauchy(nonleaf **Qq,std::pair<int, int>qSize, double **Xx,std::pair<int, int>xSize,int ifTrans,double N){
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
    nonleaf *Q = *Qq;
    double *X = *Xx;
    if(ifTrans == 0){

    }
    
    else{
        // 1st deflation permutation
        int *temp = new int[Q->TSize.first];
        for(int i = 0; i < Q->TSize.first; i++)
            temp[i] = (Q->T[i]) + 1;        
       
        arrange_elements2(&X, xSize, &(temp), Q->TSize);
        delete [] temp;      

        //double *tempX = X + (Q->n1)*xSize.second;
        double *tempXX = new double[(xSize.first - Q->n1) * xSize.second];
        memcpy(tempXX, X+(Q->n1 * xSize.second), sizeof(double)*((xSize.first - Q->n1) * xSize.second));

        bsxfun('T',&tempXX,{xSize.first - Q->n1,xSize.second},Q->v2c,Q->v2cSize);

        //eigenvalue sorting permutation
        int *tempI = new int[Q->ISize.first];        
        for(int row = 0; row < Q->ISize.first; row++)
            tempI[row] = Q->I[row] + 1;

        arrange_elements2(&tempXX, {xSize.first - Q->n1, xSize.second}, &(tempI), Q->ISize);

        delete[] tempI;
        // Givens rotation


        // second deflation permutation
        int *tempJ = new int[Q->JSize.first];
        for(int i = 0; i < Q->JSize.first; i++)
            tempJ[i] = Q->J[i] + 1;

        arrange_elements2(&tempXX, {xSize.first - Q->n1, xSize.second}, &(tempJ), Q->JSize);
        
        memcpy(X+(Q->n1 * xSize.second), tempXX, sizeof(double)*((xSize.first - Q->n1) * xSize.second));
        delete[] tempXX;
        
        int rowSize = (Q->n)-(Q->n1 + Q->n2 + 1) + 1;
        double **Qc = Q->QC;
        //tempXX = X + (Q->n1 + Q->n2 + 1) * xSize.second;
        tempXX = new double[rowSize * xSize.second];
        //for(int row = (Q->n1 + Q->n2), curr = 0; row < Q->n, curr < rowSize; row++, curr++)
        //    memcpy(tempXX + curr*(xSize.second), X+row*(xSize.second), sizeof(double)*(xSize.second));
        memcpy(tempXX, X+(Q->n1 + Q->n2)*xSize.second, sizeof(double)*(rowSize * xSize.second));
        cauchylikematvec(&Qc, (Q->qcSizes), &(tempXX), {rowSize, xSize.second}, 1);
        memcpy(X+(Q->n1 + Q->n2)*xSize.second, tempXX, sizeof(double)*(rowSize * xSize.second));
        delete[] tempXX;
       
    }
    

}


