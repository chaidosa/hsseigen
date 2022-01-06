#include<bits/stdc++.h>
#include <iostream>
#include "eigenmatrix.h"
#include "bsxfun.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}

using namespace std;

double vec_norm(vector<double> v){
    double result = 0;
    for(int i = 0; i < v.size(); ++i){
        result +=v[i]*v[i];
    }

return sqrt(result);    
}

vector<double> diff(vector<double> V){
    vector<double> result;
    for(int i = 0; i < V.size()-1; ++i){
        double temp = V[i+1]-V[i];
        result.push_back(temp);
    }
return result;
}

void rootfinder(vector<double> d,vector<double>v){
/*
%%% Input:
%%% d, v: as in secular equation
%%% N: size threshold to use fmm

%%% Ouput:
%%% x: roots of secular equation 
%%% tau: gaps
%%% org: shifts for each root
*/
    double N = 1024;
    

    int n = v.size();
    if(d.size()!=v.size()){
        cout << "not computable";
        assert(false);
    }
    //edge case
    if(n == 1){
        vector<double>x;
        vector<double>tau;
        vector<double>org;
        double temp = (v[0]*v[0]) + d[0];
        x.push_back(temp);
        temp = temp - d[0];
        tau.push_back(temp);
        org.push_back(1);
        return;
    }
   
    int C = 64;
    int record = 1;
    int alpha = 0;
    int r = 50;
    long double N0 = 1048576;
    double percent = 1;
    double FMM_ITER = ceil(log(n)/log(2))-6 > 5?ceil(log(n)/log(2))-6:5;
    int Flag = 1;
    int MAX_ITER = 100;
   
    if(record){
        vector<double> iter(n);
        vector<double> residual(n);
        vector<double> erretms(n);
    }
    double v_norm = vec_norm(v);
    double rho    = 1 / (v_norm*v_norm);

    for(int i = 0; i < v.size(); ++i){
        v[i] = v[i] / v_norm;
    }
    vector<double> d0 = diff(d);

    vector<double> x;
    for(int i = 0; i < n-1; ++i){
        double temp = (d[i] + d[i+1]) / 2;
        x.push_back(temp);
    }

    vector<double> org;
    for(int i=0; i < n-1; ++i){
        org.push_back(i);
    }

   double *f0;
   int f0_size;
   double *v2_arr;
   v2_arr = new double[v.size()];
   for(int i = 0; i < v.size(); ++i){
        v2_arr[i] = v[i]*v[i];
    }

   if(n >= N){
        //[fl, fu, nflops1] = trifmm1d_local_shift(r, x, d, v.^2, d0/2, org, 1);
        //f0 = rho - fl - fu;
   }
   else{
       //  K = d(org) - d.';
       //  K = bsxfun(@plus, K, d0/2);
       //  K = 1 ./ K;
        int kRows = org.size();
        int kCols  = d.size(); 
        std::pair<int,int> Ksize = {kRows, kCols};
        double *K = new double[kRows*kCols];
        for(int row = 0; row < kRows; ++row){
            for(int col=0; col < kCols; ++col){
                K[col+row*(kCols)] = d[row] - d[col];
            }
        }
        double *tempD = new double[d0.size()];
        for(int i = 0; i < d0.size(); ++i){
            tempD[i] = d[i] / 2;
        }

        bsxfun('P',&K,{kRows,kCols},tempD,{d0.size(),1});
        
        for(int row_col = 0; row_col < (kRows*kCols); ++row_col){
            K[row_col] = 1 / K[row_col];
        }        
        // f0 = rho - K * v.^2;
        f0 = new double[kRows]; 
        f0_size = kRows;
        memcpy(f0, 0, sizeof(double)*f0_size);   
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1.0, K, kCols, v2_arr, 1, 0.0, f0, 1);

        for(int i = 0; i < kRows; ++i){
            f0[i] = rho - f0[i];
        }
    }

    // h = 2 * diff(v.^2) ./ d0;
    vector<double>tempvSqr(v2_arr,v2_arr+v.size());
    vector<double> h = diff(tempvSqr);
    for(int i = 0; i < h.size(); ++i){
        h[i] = 2*h[i] / d0[i];
    }
    //g = f0 - h;
    vector<double>g;
    for(int i = 0; i < f0_size; ++i){
        g.push_back((f0[i]-h[i]));
    }
    // org = [1:n-1]' + (f0<0);
    for(int i = 0; i < org.size(); ++i){
        if(f0[i] < 0){
            org[i] = org[i] + 1;
        }
    }    

    //Upper bound and lower bound
    vector<double> d1(d.begin(),d.end()-1);
    vector<double> d2(d.begin()+1,d.end());
    vector<double> v1(v.begin(),v.end()-1);
    vector<double> v2(v.begin()+1, v.end());
    vector<double> xub(n-1,0);
    vector<double> xlb(n-1,0);
    vector<int> I1;
    vector<int> I2;
    for(int i = 0; i < f0_size; i++){
        if(f0[i]>=0)
            I1.push_back(i);
        else
            I2.push_back(i);
    }

    for(int i = 0; i < I1.size(); ++i){
        xub[I1[i]] = (d2[I1[i]] - d1[I1[i]]) / 2;  
    }
    for(int i = 0; i < I2.size(); ++i){
        xlb[I2[i]] = (d1[I2[i]] - d2[I2[i]]) / 2;  
    }
    vector<double>xlb0(xlb.begin(),xlb.end());
    vector<double>xub0(xub.begin(),xub.end());

    //inital guess
    //a = (2*(f0>=0) - 1) .* d0 .* g + v(1:n-1).^2 + v(2:n).^2;    
    vector<double>a;
    for(int i = 0; i < n-1; ++i){
        double temp_first_part = 0;
        double temp_second_part = 0;
        if(f0[i]>=0){
            temp_first_part = 1;
        }
        else{
            temp_first_part = -1;
        }
        temp_first_part+=d0[i]*g[i];
        temp_second_part = tempvSqr[i]+tempvSqr[i+1];
        temp_first_part+= temp_second_part;
        a.push_back(temp_first_part);
    }
    //b = ((f0>=0) .* v(1:n-1).^2 - (f0<0) .* v(2:n).^2) .* d0;
    vector<double>b;
    for(int i = 0; i < n-1; ++i){
        double temp_first_part = 0;
        double temp_second_part = 0;        
        if(f0[i] >= 0){
            temp_first_part = tempvSqr[i];
            temp_second_part = 0;

        }
        else{
            temp_first_part = 0;
            temp_second_part = tempvSqr[i+1]*d0[i];
        }
        temp_first_part = temp_first_part - temp_second_part;
        b.push_back(temp_first_part);
    }

    vector<double>tau(n-1);    
    std::fill(tau.begin(), tau.end(), 0);
    I1.resize(0);
    I1.shrink_to_fit();
    I2.resize(0);
    I2.shrink_to_fit();

    for(int i = 0; i < a.size(); ++i){
        if(a[i]>0){
            I1.push_back(i);
        }
        else{
            I2.push_back(i);
        }
    }
    //tau(I1) = 2*b(I1) ./ (a(I1) + sqrt(abs(a(I1).^2 - 4*b(I1).*g(I1))));    
    for(int i = 0; i < I1.size(); ++i){
        tau[I1[i]] = 2*b[I1[i]] / (a[I1[i]] + std::sqrt(abs( a[I1[i]]*a[I1[i]] - 4*b[I1[i]] * g[I1[i]])));
    }
    //tau(I2) = (a(I2) - sqrt(abs(a(I2).^2 - 4*b(I2).*g(I2)))) ./ (2*g(I2));
    for(int i = 0; i < I2.size(); ++i){
        tau[I2[i]] = (a[I2[i]] - std::sqrt(abs( a[I2[i]]* a[I2[i]] - (4 * b[I2[i]]*g[I2[i]])))) / (2* g[I2[i]]);
    }

    //bound check

}
