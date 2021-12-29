#include "cauchylikematvec.h"
#include "bsxfun.h"
#include <string.h>
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}
void cauchylikematvec(double **Qc, std::pair<int,int>*qcSizes,double *X, std::pair<int,int>xSize,int ifTrans, double N=1024){
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
    double alpha = 1.0;
    double beta  = 0.0;

    double *v = Qc[0];
    double *s = Qc[1];
    double *d = Qc[2];
    double *lam = Qc[3];
    double *tau = Qc[4];
    double *org = Qc[5];
    int n = qcSizes[0].first;
    int r =50;

    if(n < N){
        if(ifTrans == 0){

        }
        else{
            double *d_org = new double[qcSizes[2].first];
            arrange_elements(d,qcSizes[2],org,qcSizes[5],d_org);
            double *S = new double[(qcSizes[2].first)*(qcSizes[2].first)];  

            //d_org - d.'
            for(int row = 0;row < qcSizes[2].first; row++){
                for(int col=0; col < qcSizes[2].first; col++){
                    S[col + row*(qcSizes[2].first)] = d_org[row] - d[col];
                }
            }

            bsxfun('P',S,{qcSizes[2].first,qcSizes[2].first},tau,qcSizes[5]);
            //S = -1 ./S;
            for(int row = 0;row < qcSizes[2].first; row++){
                for(int col=0; col < qcSizes[2].first; col++){
                    S[col + row*(qcSizes[2].first)] = -1 / S[col + row*(qcSizes[2].first)];
                }
            }
            bsxfun('T',S,{qcSizes[2].first,qcSizes[2].first},s,qcSizes[1]);
            bsxfun('T',S,{qcSizes[2].first,qcSizes[2].first},v,{qcSizes[0].second,qcSizes[0].first}); 
            //Y = S*X
            double *Y = new double[qcSizes[2].first*xSize.second];
            memset(Y,0,sizeof(double)*(qcSizes[2].first*xSize.second));
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,qcSizes[2].first,xSize.second,qcSizes[2].second,alpha,S,qcSizes[2].second,X,xSize.second,beta,Y,xSize.second);
            delete [] X;
            X = Y;
            delete [] d_org;
            delete [] S;
        }
    }
    else{
        //FMM
    }

}