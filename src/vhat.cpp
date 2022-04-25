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

std::vector<double> vhat(std::vector<double>& d, std::vector<double>& lam, const int *org, int org_size, std::vector<double>& tau, std::vector<double>& w, double N){

    int n = lam.size();
    int r = 50;
    double *Dlam,*D;
    std::vector<double>v(n);

    if(n < N){
        int dRows = d.size();
        int dCols = org_size;
        Dlam = new double[dRows*dCols];
        
        for(int row = 0; row < dRows; row++)
            for(int col = 0; col < dCols; col++)
                Dlam[col + row*dCols] = d[row] - d[org[col]];

                
        //bsxfun('M','R', &Dlam, {dRows, dCols}, &tau[0], {1, tau.size()});
        for(int row = 0 ; row < dRows; row++){
                for(int col = 0; col<dCols; col++){
                    Dlam[col+row*dCols] = Dlam[col+row*dCols]-tau[col];
                }
            }


        // D = d - d.' and D(1:(n+1):end) = ones(n,1)
        D = new double[dRows*dRows];
        for(int row = 0; row < dRows; row++){
            for(int col = 0; col < dRows; col++){
                D[col + row*dRows] = d[row] - d[col];
                if(row == col)
                    D[col + row*dRows] = 1;
            }
        }


        // v = (log(abs(Dlam))-log(abs(D))) * ones(n,1);
        for(int row = 0; row < dRows; row++){
            for(int col = 0; col < dCols; col++){
                v[row] += std::log(std::abs(Dlam[col + row*dCols])) -std::log(std::abs(D[col + row*dCols]));
            }
        }


        //v = exp(1/2 * v); and v(w < 0) = -v(w < 0)
        for(int row = 0; row < v.size(); row++){
            double temp = v[row]/2;
            v[row] = std::exp(temp);
        }     

    }

    else{
        //fmm1d_local_shift2()
        //fmm1d_local_shift()
    }

    delete [] Dlam;
    delete [] D;     
   
    return v;
}
