#include<bits/stdc++.h>
#include<fstream>
#include<iomanip>
#include<algorithm>
#include <iostream>
#include<string.h>
#include "bsxfun.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}

void colnorms(std::vector<double>& d, std::vector<double>& lam, std::vector<double>& tau, std::vector<double>& org, std::vector<double>& v, double N = 1024){
/*
%%% Input:
%%% d, lam, tau, org, v : from secular.m and rootfinder.m
%%% N: size threshold to use fmm

%%% Output
%%% s: norm of each columns of ( v_i / (d_i - lam_j) )
*/ 
    int n = lam.size();
    int r = 50; 
    double *S, *s;

    if(n < N){
        int sRows = org.size();
        int sCols = d.size();
        S = new double[sRows*sCols];

        for(int row = 0; row < sRows; row++){
            for(int col = 0; col < sCols; col++){
                S[col + row *sCols] = d[(int)org[row]] - d[col];
            }
        }

        bsxfun('P',&S,{sRows, sCols}, &tau[0], {tau.size(), 1});

        for(int row_col = 0; row_col < (sRows*sCols); row_col++)
            S[row_col] = 1 / S[row_col];

        s = new double[sRows];
        memset(s, 0 ,sizeof(double)*sRows*sCols);

        double *temp_vSqr = new double[v.size()];
        for(int i = 0; i < v.size(); ++i)
            temp_vSqr[i] = v[i]*v[i];
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sRows, 1, sCols, 1.0, S, sCols, temp_vSqr, 1, 0.0, s, 1);

        for(int i = 0; i < sRows; i++)
            s[i] = 1 / std::sqrt(s[i]);

    }
    else{
        //FMM 1D local shift
    }

}