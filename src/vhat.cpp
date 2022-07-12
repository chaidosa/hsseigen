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
#include "fmm1d_local_shift_2.h"
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
        assert(false);
	assert(d.size() == n);
	assert(w.size() == n);
        
	double* ptrD=d.data();
	double* gap=tau.data();
	std::vector<double> e1(n,1);
	std::vector<double> zeros(n,0);
	//Vec.
	double* v0 = fmm1d_local_shift_2(r, ptrD, lam, e1.data(), gap, org, 3, d.size(), n);
        double* vd = fmm1d_local_shift(r, ptrD, ptrD, e1.data(), zeros.data(), org, 3, d.size(), d.size());
	for(int i=0;i<n;i++)
		v[i]=pow(0.5, v0[i]-vd[i]);
 
	//Vec:
	for(int i=0;i<n;i++)
		if(w[i] <0)
			v[i] = -1 * v[i];
    	delete [] v0;
	delete [] vd;
    }
    return v;
}
