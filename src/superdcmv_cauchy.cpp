#include "superdcmv_cauchy.h"
#include <iostream>
#include <utility>
#include <string.h>
#include "eigenmatrix.h"
#include "bsxfun.h"
#include "cauchylikematvec.h"
#include <fstream>
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

#include <iomanip>

void printX(string file, double *X){
    std::ofstream OutFile(file);
    for(int i=0; i<32; i++){
        OutFile << setprecision(20) << X[i] << "\n";
    }
    OutFile.close();
}

void loadX(string file, double *X){
    std::ifstream InFile(file);
    for (int i = 0; i < 32; i++)
    {
        InFile >> setprecision(20) >> X[i];
    }
    InFile.close();
    
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
	//orthogonal Cauchy eignematrix
        int load=0;
        if(load){loadX("X_cauchy_cauchylikein.txt", X);}
	int rowSize = (Q->n)-(Q->n1 + Q->n2 + 1) + 1;
        double **Qc = Q->QC;        
        double* ret = cauchylikematvec(Qc, (Q->qcSizes), Q->Org, X+(Q->n1 + Q->n2)*xSize.second, {rowSize, xSize.second}, ifTrans, N);
        if(load){printX("X_cauchy_cauchylike_T=0.txt", ret);}
        memcpy(X+(Q->n1 + Q->n2)*xSize.second, ret, sizeof(double)*(rowSize * xSize.second));
        load = 0;


	delete[] ret;

	//2nd deflation permutation
        double *tempXX = new double[(xSize.first - Q->n1) * xSize.second];
        memcpy(tempXX, X+(Q->n1 * xSize.second), sizeof(double)*((xSize.first - Q->n1) * xSize.second));
	int *tempJ = new int[Q->JSize.first];
        for(int i = 0; i < Q->JSize.first; i++)
            tempJ[i] = Q->J[i] + 1;
	/*calls _lapmr, which rearranges the rows of the m-by-n matrix tempXX as specified by the permutation tempJ(1),tempJ(2),...,tempJ(m) of the integers 1,...,m.. If the last argument is true (by default) then tempXX[i]=tempXX[J[i]] otherwise tempXX[J[i]]=tempXX[i]*/
	assert(Q->JSize.first == (xSize.first - Q->n1));
        arrange_elements2(tempXX, {xSize.first - Q->n1, xSize.second}, (tempJ), Q->JSize, false);

        // Givens rotation. 
	/*Optimizing cblas_dgemm usage without using temporary array tempArr. Using ldb argument of cblas_dgemm and inspecting p and j appropriately. e.g. if p<j, then ldb=(j-p)*xSize.second and tempXX+p*xSize.second is Matrix B. If j<p, then matrix B=tempXX+j*xSize.second and ldb=(p-j)*xSize.second  */

        //double *tempArr = new double[2*xSize.second];
        double *MulR = new double[2*xSize.second];
	assert(((Q->GSize.first) % 4) == 0);
        for(int l = (Q->GSize.first)-1; l >=0 ; l = l - 4){
            int p =(int)Q->G[l-3];
            int j =(int)Q->G[l-2];
            double c = Q->G[l-1];
            double s = Q->G[l];            
            
	    assert(j == p+1); //delete this after a single round of testing.
	    int ldb=(j-p)*xSize.second;
	    double* matrixB=tempXX+p*xSize.second;
	    if(j<p) {
	    	ldb=(p-j)*xSize.second;
	    	matrixB=tempXX+j*xSize.second;
	    }

            //memcpy(tempArr, tempXX + p*xSize.second, sizeof(double)*(xSize.second));
            //memcpy(tempArr+xSize.second, tempXX + j*xSize.second, sizeof(double)*(xSize.second));
            double GivensArr[2*2] = {c, s, -s, c};
           
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, xSize.second, 2, 1, GivensArr, 2, matrixB, ldb, 0.0, MulR, xSize.second);
            
            memcpy(tempXX + p *xSize.second, MulR, sizeof(double)*(xSize.second));
            memcpy(tempXX + j *xSize.second, MulR + xSize.second, sizeof(double)*(xSize.second));
        }

        //delete [] tempArr;
        delete [] MulR;
 	
	//eigenvalue sorting permutation
	assert(Q->ISize.first == (xSize.first - Q->n1));
	//if this assertion holds, then TODO: remove the temporary memory creation for tempI and reuse tempJ created earlier.
        int *tempI = new int[Q->ISize.first];        
        for(int row = 0; row < Q->ISize.first; row++)
            tempI[row] = Q->I[row] + 1;

        arrange_elements2(tempXX, {xSize.first - Q->n1, xSize.second}, (tempI), Q->ISize, false);
        delete[] tempI;



	//conjugate normalizer. TODO: when v2c is complex, take conjugate of v2c and pass as argument.
        bsxfun('T',&tempXX,{xSize.first - Q->n1,xSize.second},Q->v2c,Q->v2cSize);

	// 1st deflation permutation
        memcpy(X+(Q->n1 * xSize.second), tempXX, sizeof(double)*((xSize.first - Q->n1) * xSize.second));
        delete[] tempXX;

	assert(Q->TSize.first == xSize.first);
        int *temp = new int[Q->TSize.first];
        for(int i = 0; i < Q->TSize.first; i++)
            temp[i] = (Q->T[i]) + 1;        
       
        arrange_elements2(X, xSize, (temp), Q->TSize, false);
        delete [] temp;      
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
           
           assert(j == p+1);
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
        
        int print = 0;
        if (print){printX("X_cauchy_cauchylikein.txt", X);}
        int rowSize = (Q->n)-(Q->n1 + Q->n2 + 1) + 1;
        double **Qc = Q->QC;        
        double* ret = cauchylikematvec(Qc, (Q->qcSizes), Q->Org, X+(Q->n1 + Q->n2)*xSize.second, {rowSize, xSize.second}, 1, N);
        if (print){printX("X_cauchy_cauchylike_T=1.txt", ret);}
        print = 0;
        memcpy(X+(Q->n1 + Q->n2)*xSize.second, ret, sizeof(double)*(rowSize * xSize.second));
	delete[] ret;
    }
    
return X;
}


// int rowSize = (Q->n)-(Q->n1 + Q->n2 + 1) + 1;
//         double **Qc = Q->QC;        
//         double* ret = cauchylikematvec(Qc, (Q->qcSizes), Q->Org, X+(Q->n1 + Q->n2)*xSize.second, {rowSize, xSize.second}, ifTrans, N);
//         // int print = 0;
//         memcpy(X+(Q->n1 + Q->n2)*xSize.second, ret, sizeof(double)*(rowSize * xSize.second));