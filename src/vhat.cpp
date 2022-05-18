/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% compute vhat by Lowner's formula %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include<bits/stdc++.h>
#include<fstream>
#include<iomanip>
#include <iostream>
#include<fstream>
#include<iomanip>
#include "bsxfun.h"
#include "vhat.h"

double * vhat(std::vector<double>& d, double* lam, const int *org, int org_size, std::vector<double>& tau, std::vector<double>& w, double N){

    int n = org_size;
    int r = 50;   
    //std::vector<double>v(n);
    double *v = new double[n];
    if(n < N){
        int dRows = d.size();
        int dCols = org_size;
        for(int row = 0; row < dRows; row++){
            v[row] = 0;
            for(int col = 0; col < dCols; col++){
                double temp = d[row] - d[org[col]];
                temp = temp - tau[col];
                //Dlam[col + row*dCols]
                //D[col + row*dRows] = d[row] - d[col];
                double temp2 = d[row] - d[col];
                if(row == col)
                    temp2 = 1;

                v[row] += std::log(std::abs(temp)) -std::log(std::abs(temp2));

            }

            double temp3 = v[row]/2;
            v[row] = std::exp(temp3);
        }
        for(int itr = 0; itr < w.size(); itr++){
            if(w[itr] < 0)
                v[itr] = -v[itr];
        }              
    }
    else{
        //fmm1d_local_shift2()
        //fmm1d_local_shift()
    }
 
   
    return v;
}
