#include<bits/stdc++.h>
#include<fstream>
#include<iomanip>
#include<algorithm>
#include <iostream>
#include<string.h>
#include "bsxfun.h"
#include "colnorms.h"
#include "fmm1d_local_shift_2.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}

double* colnorms(std::vector<double>& d, double* lam, std::vector<double>& tau, const int *org, int org_size, double *v, double N){
/*
%%% Input:
%%% d, lam, tau, org, v : from secular.m and rootfinder.m
%%% N: size threshold to use fmm

%%% Output
%%% s: norm of each columns of ( v_i / (d_i - lam_j) )
*/ 
    int n = org_size;
    int r = 50; 
    double *S=NULL, *s=NULL;
    int sSize=0;

    double *temp_vSqr = new double[org_size];
    for(int i = 0; i < org_size; ++i)        
	temp_vSqr[i] = v[i]*v[i];
    
    //if(n < N)
    {
        int sRows = org_size;
        int sCols = d.size();
        S = new double[sRows*sCols];

        for(int row = 0; row < sRows; row++){
            for(int col = 0; col < sCols; col++){
                double temp = d[org[row]] - d[col] + tau[row];
                S[col + row *sCols] = 1 / (temp * temp);
            }
        }
        
        //bsxfun('P',&S,{sRows, sCols}, &tau[0], {tau.size(), 1});        

        //for(int row_col = 0; row_col < (sRows*sCols); row_col++)
        //    S[row_col] = 1 / (S[row_col]*S[row_col]);

        
        s = new double[sRows];
        sSize = sRows;
        memset(s, 0 ,sizeof(double)*sRows);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sRows, 1, sCols, 1.0, S, sCols, temp_vSqr, 1, 0.0, s, 1);

        for(int i = 0; i < sRows; i++)
            s[i] = 1 / std::sqrt(s[i]);

        delete [] S;
    }
#if 0
    else{
        //FMM 1D local shift
    	s=fmm1d_local_shift(r,lam,d.data(),temp_vSqr,tau.data(), org, 2, org_size, org_size);
    	for(int i=0;i<org_size;i++){
		assert(s[i] >=0);
	    s[i]=1/std::sqrt(s[i]);
    	}
    }
#endif
    delete [] temp_vSqr;

    //std::vector<double>result(s, s+sSize);
    return s;
}
