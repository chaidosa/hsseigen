#include "superdcmv_cauchy.h"
#include <iostream>
#include <utility>
#include <string.h>
#include "eigenmatrix.h"
#include "bsxfun.h"
#include "cauchylikematvec.h"
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
        double *temp;
        arrange_elements2(&X,xSize,&(Q->T),Q->TSize,&temp);

        delete [] X;
        X = temp;
        temp = NULL;

        xSize = {Q->TSize.first,xSize.second};

        //conjugate normalizermalloc(): corrupted top size

        double *tempX = X + (Q->n1)*xSize.second;
        bsxfun('T',&tempX,{xSize.first - Q->n1,xSize.second},Q->v2c,Q->v2cSize);

        //eigenvalue sorting permutation
        double *tempI = new double[Q->ISize.first];
        for(int row = 0; row < Q->ISize.first; row++)
            tempI[row] = Q->I[row] + Q->n1;

        double *result_eigen_permutation;

        arrange_elements2(&tempX,{xSize.first - Q->n1,xSize.second},&tempI,{Q->ISize.first, Q->ISize.second},&result_eigen_permutation);
        memcpy(tempX,result_eigen_permutation,sizeof(double)*((xSize.first - Q->n1)*xSize.second));
        delete [] tempI;
        //Givens rotation




        //2nd deflation permutation
        double *res_2_deflation;

        double *tempJ = new double[Q->JSize.first];
        for(int row = 0; row < Q->JSize.first; row++)
            tempJ[row] = Q->J[row] + Q->n1;

        arrange_elements2(&tempX,{xSize.first - Q->n1,xSize.second},&tempJ,{Q->JSize.first, Q->JSize.second},&res_2_deflation);
        delete [] tempJ;
         
        memcpy(tempX,res_2_deflation,sizeof(double)*((xSize.first - Q->n1)*xSize.second));
        
        double **Qc = Q->QC;

        int tempRows = Q->n1+Q->n2;
        int tempRowe = Q->n;
        int tempcols = 0;
        int tempcole = xSize.second;
        double *tempp = new double[(tempRowe-tempRows)*(tempcole)];
        memcpy(tempp, X+tempRows, sizeof(double)*((tempRowe-tempRows)*(tempcole)));
        cauchylikematvec(&Qc, Q->qcSizes, &tempp, {(tempRowe-tempRows), tempcole}, 1);        
        memcpy(X+tempRows, tempp, sizeof(double)*((tempRowe-tempRows)*(tempcole)));
        std::cout<<"Success\n";
        *Xx = X;
        delete [] tempp;
        delete [] res_2_deflation;
        delete [] result_eigen_permutation;
    }
    

}


