#include "cauchylikematvec.h"
#include "bsxfun.h"
#include <string.h>
#include "fmm1d_local_shift_2.h"
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
double *cauchylikematvec(double **Qcd, std::pair<int,int>*qcSizes, const int *org, double *Xx, std::pair<int,int>xSize,int ifTrans, double N){
/*
    %%% Input:
    %%% Qc: {v,s,d,lam} stores the cauchylike matrix  v_i * s_j / (d_i - lam_j) 
    %%% X: vectors to be multiplied
    %%% ifTrans = 0: not transpose
    %%% ifTrans = 1: transpose
    %%% N: size threshold to use fmm

    %%% Output
    %%% t = 0: Y = Q * X 
    %%% t = 1: Y = Q^T * X;
*/
    double *Y=NULL;//return value

    double **Qc = Qcd;
    double *X    = Xx;
    double alpha = 1.0;
    double beta  = 0.0;

    double *v    = Qc[0];
    double *s    = Qc[1];
    double *d    = Qc[2];
    double *lam  = Qc[3];
    double *tau  = Qc[4];
    //double *org  = Qc[5]; //here change to be made
    //int *orgd; // change the above org
    int n        = qcSizes[0].first;
    int r        = 50;

    if(n < N)
    {
        if(ifTrans == 0){
            assert(false);
        }

        else
        {
            double *S = new double[(qcSizes[5].first)*(qcSizes[2].first)];  
            int sRows = qcSizes[2].first;
            int sCols = qcSizes[5].first;
            //d_org - d.'
            for(int row = 0;row < sRows; row++){
                for(int col = 0; col < sCols; col++){
                    double temp = d[org[row]] - d[col];
                    //S[col + row*(qcSizes[2].first)] = d[org[row]] - d[col];
                    temp += tau[row];
                    temp = -1 / temp;
                    temp = temp * s[row];
                    temp = temp * v[col];
                    S[col + row*(qcSizes[2].first)] = temp;
                }
            }
       
            Y = new double[qcSizes[5].first*xSize.second];
            //memset(Y,0,sizeof(double)*(qcSizes[5].first*xSize.second));
            
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,qcSizes[5].first,xSize.second,qcSizes[2].first,alpha,S,qcSizes[2].first,X,xSize.second,beta,Y,xSize.second);           
            delete [] S; 
        }
    }

    else{
        //FMM
        if(ifTrans == 0){
            bsxfun('T', &X, xSize, s, qcSizes[1]);
            Y = fmm1d_local_shift_2(r, d, lam, X, tau, org, 1, qcSizes[2].first, qcSizes[3].first, 1);
            bsxfun('T', &Y, qcSizes[2], v, qcSizes[0]);
        }
        else{
            bsxfun('T', &X, xSize, v, qcSizes[0]);
            Y = fmm1d_local_shift(r, lam, d, X, tau, org, 1, qcSizes[3].first, qcSizes[2].first, 1);
	    
	    //Vec:
	    for(int i=0;i<qcSizes[3].first;i++)
		    Y[i] = -1 * Y[i];

            bsxfun('T', &Y, qcSizes[3], s, qcSizes[1]);
        }
    }
    return Y;
}
