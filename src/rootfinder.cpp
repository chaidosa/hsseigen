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
    double alpha = 0;
    int r = 50;
    long double N0 = 1048576;
    double percent = 1;
    double FMM_ITER = ceil(log(n)/log(2))-6 > 5?ceil(log(n)/log(2))-6:5;
    int Flag = 1;
    int MAX_ITER = 100;
    double eps = 2.2204e-16;
   // if(record){
        vector<double> iter(n);
        vector<double> residual(n);
        vector<double> erretms(n);
  //  }
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
   int f0_size = 0;
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
        
        double *K = new double[kRows*kCols];

        for(int row = 0; row < kRows; ++row){
            for(int col=0; col < kCols; ++col){
                K[col+row*(kCols)] = d[row] - d[col];
            }
        }

        double *tempD = new double[d0.size()];

        for(int i = 0; i < d0.size(); ++i){
            tempD[i] = d0[i] / 2;
        }

        bsxfun('P',&K,{kRows,kCols},tempD,{d0.size(),1});
        
        for(int row_col = 0; row_col < (kRows*kCols); ++row_col){
            K[row_col] = 1 / K[row_col];
        }

        // f0 = rho - K * v.^2;
        f0 = new double[kRows]; 
        f0_size = kRows;
        memset(f0, 0, sizeof(double)*kRows); 

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1.0, K, kCols, v2_arr, 1, 0.0, f0, 1);

        for(int i = 0; i < kRows; ++i){
            f0[i] = rho - f0[i];
        }

        delete[] K;
        K = NULL;
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

        temp_first_part *= d0[i]*g[i];
        temp_second_part = tempvSqr[i]+tempvSqr[i+1];
        temp_first_part += temp_second_part;
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
            temp_second_part = tempvSqr[i+1];

        }

        temp_first_part = (temp_first_part - temp_second_part)*d0[i];
        b.push_back(temp_first_part);
    }

    vector<double>tau(n-1);    
    
    I1.resize(0);
    I1.shrink_to_fit();
    I2.resize(0);
    I2.shrink_to_fit();

    for(int i = 0; i < a.size(); ++i){

        if(a[i]>0)
            I1.push_back(i);        

        else
            I2.push_back(i);
        
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
    std::vector<int>I_B;
    for(int i = 0; i < xlb.size(); ++i){

        if(tau[i]<=xlb[i] || tau[i]>=xub[i])
            I_B.push_back(i);

    }

    if(I_B.size() != 0){

        for(int i = 0; i < I_B.size(); ++i){
            tau[I_B[i]] = (xub[I_B[i]] + xlb[I_B[i]])/2;
        }

    }
    
    for(int i = 0; i < tau.size(); ++i){
        x[i] = tau[i] + d[(int)org[i]];
    }

// iterate n-1 roots simultaneously
    int iter_ct = 0; //iteration count
    vector<double> swtch(n-1,1);
    vector<double> pref;
    vector<double> f_outer;
    std::vector<double> J2;
    double *f, *psi, *phi, *dpsi, *dphi, *df;     
        int pd_size = 0;
    while (iter_ct < FMM_ITER && Flag){
    
        if(iter_ct > 0)
            pref = f_outer;
    
        //bound check with bisection safeguard
        I_B.clear();
        I_B.resize(0);
        I_B.shrink_to_fit();

        for(int i = 0; i < tau.size(); i++){
            if(tau[i] == 0 || (tau[i]*xub[i]) < 0 || (tau[i]*xlb[i])<0 || isnan(tau[i]) || tau[i] < xlb0[i] || tau[i] > xub0[i])
                I_B.push_back(i);
        }

        if(I_B.size() != 0){
            for(int i = 0; i < I_B.size(); ++i){
                tau[I_B[i]] = (xub[I_B[i]] + xlb[I_B[i]])/2;
                x[I_B[i]] = tau[I_B[i]] + d[(int)org[I_B[i]]];
            }
        }

        //secular function evaluation
        
        
        if(n >= N){
            //fmm
        }

        else{

            int kRows = org.size();
            int kCols  = d.size();             
            double *K = new double[kRows*kCols];

            for(int row = 0; row < kRows; ++row){
                for(int col=0; col < kCols; ++col){
                    K[col+row*(kCols)] = d[row] - d[col];
                }
            }

            double* temp_tau = &tau[0];

            //double* temp_org;
            //memcpy(temp_org,temp_o,sizeof(int)*org.size());
            bsxfun('P',&K,{kRows,kCols},temp_tau,{tau.size(),1});

            for(int i = 0; i < (kRows*kCols); ++i)
                K[i] = 1/K[i];

            double *K1 = new double[kRows*kCols];
            double *K2 = new double[kRows*kCols];

            memcpy(K1,K,sizeof(double)*kRows*kCols);
            memcpy(K2,K,sizeof(double)*kRows*kCols);

            //K1 = tril(K), K2 = triu(K,1)
            for(int row = 0; row < kRows; row++){
                for(int col = 0; col < kCols; col++){
                    if(row < col)
                        K1[col + row*(kCols)] = 0;
                    else
                        K2[col + row*(kCols)] = 0;
                }
            }

            psi = new double[kRows];
            phi = new double[kRows];

            memset(psi,0,sizeof(double)*kRows);
            memset(phi,0,sizeof(double)*kRows);

            double *temp = new double[v.size()];
            for(int i = 0; i < v.size(); ++i)
                temp[i] = v[i]*v[i];

            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,kRows,1,kCols,1,K1,kCols,temp,1,0,psi,1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,kRows,1,kCols,1,K2,kCols,temp,1,0,phi,1);
            
            f = new double[kRows];

            for(int i = 0; i < kRows; i++)
                f[i] = rho-psi[i]-phi[i];
            
            double *dK1 = K1;
            double *dK2 = K2;

            //no of flops for multiplication can be reduced here.
            for(int row = 0; row < kRows; row++){
                for(int col = 0; col < kCols; col++){

                    dK2[col + row*kCols] = dK2[col + row*kCols]*dK2[col + row*kCols];        
                    dK1[col + row*kCols] = dK1[col + row*kCols]*dK1[col + row*kCols];

                }
            }

            dpsi = new double[kRows];
            dphi = new double[kRows];

            memset(dpsi,0,sizeof(double)*kRows);
            memset(dphi,0,sizeof(double)*kRows);

            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,kRows,1,kCols,1,dK1,kCols,temp,1,0,dpsi,1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,kRows,1,kCols,1,dK2,kCols,temp,1,0,dphi,1);

            df = new double[kRows];
            for(int i = 0; i < kRows; ++i){
                df[i] = dpsi[i] + dphi[i];
            }

            pd_size = kRows;

            delete[] dK1;
            delete[] dK2;
            delete[] K;

        // end of else statement
        }

        //adaptive FMM iterations        
        double *res = new double[pd_size];

        for(int i = 0 ; i < pd_size; ++i)
            res[i] = rho + std::abs(psi[i]) + std::abs(phi[i]);

        double nJ2 = 0;
        if(iter_ct > 0)
            nJ2 = J2.size();        
        
        for(int i = 0; i < pd_size; ++i){
            double max = std::max(res[i],alpha);
            max = C*n*eps*max;

            if(std::abs(f[i]) > max)
                J2.push_back(i);         
        
        }
        
        if(iter_ct > 0 && (J2.size() <= 0.01*n || std::abs(nJ2 - J2.size()) < 0.01 *n))
            Flag = Flag - 1;

        if(iter_ct == 5 || iter_ct < 5 && Flag == 0)
            percent = J2.size() / (n-1);
        
        //update  upper and lower bounds
        I_B.clear();
        I_B.resize(0);
        I_B.shrink_to_fit();        
        vector<int> II;
        for(int i = 0; i < pd_size; ++i){            
            if(f[i] < 0)
                I_B.push_back(i);
            else
                II.push_back(i);
        }

        for(int i = 0; i < I_B.size(); ++i)
            xlb[I_B[i]] = std::max(tau[I_B[i]],xlb[I_B[i]]);

        for(int i = 0; i < II.size(); ++i)
            xub[II[i]] = std::min(tau[II[i]],xub[II[i]]);

        //swtch = 1 : fixed weight method
        //swtch = 0 : middle way method
        vector<int> Iswtch;
        if(iter_ct > 0){
            for(int i = 0; i < pd_size; ++i){
                if(pref[i]*f[i] > 0 & (std::abs(pref[i]) > abs(f[i]) / 10))
                    Iswtch.push_back(i);
            }
            for(int i = 0; i < Iswtch.size(); ++i){
                swtch[Iswtch[i]] = -swtch[Iswtch[i]] + 1;
            }
        }

        //quadratic equation coefficients
        a.clear();
        a.resize(0);
        a.shrink_to_fit();
        b.clear();
        b.resize(0);
        b.shrink_to_fit();
        for(int i = 0; i < pd_size; ++i){
            double temp_a = ( (d[i]-x[i]) + (d2[i]-x[i]) ) * f[i] - (d1[i]-x[i]) * (d2[i]-x[i]) * df[i];
            a.push_back(temp_a);

            double temp_b = (d1[i]-x[i]) * (d2[i]-x[i]) * f[i];
            b.push_back(temp_b);
        }
        vector<double>c(n-1);
        if(iter_ct > 0){
            vector<double> c0;
            vector<double> c1;
            for(int i = 0; i < pd_size; ++i){
                double temp = -(d1[i]-x[i]) * dpsi[i] - (d2[i]-x[i]) * dphi[i] + f[i];
                c0.push_back(temp);
                temp = 0;
                if(f0[i]>=0)
                    temp = -1*( (d2[i] - x[i]) * df[i] + (v1[i]*v1[i]) * (d1[i]- d2[i]) / ((d1[i]-x[i])*(d1[i]-x[i])));
                else
                    temp = -1*((d1[i] - x[i]) * df[i] + (v2[i]*v2[i])*(d2[i]-d1[i]) / ((d2[i]-x[i])*(d2[i]-x[i])));

                temp = temp + f[i];
                c1.push_back(temp);
            }
            
            for(int i = 0; i < swtch.size(); ++i){
                if(swtch[i] == 0)
                    c[i] = c0[i];
                else
                    c[i] = c1[i];

            }
        }

        else{
            for(int i = 0; i < pd_size; ++i){
                double temp;
                if(f0[i]>=0)
                    temp = -1*( (d2[i] - x[i]) * df[i] + (v1[i]*v1[i]) * (d1[i]- d2[i]) / ((d1[i]-x[i])*(d1[i]-x[i])));
                else
                    temp = -1*((d1[i] - x[i]) * df[i] + (v2[i]*v2[i])*(d2[i]-d1[i]) / ((d2[i]-x[i])*(d2[i]-x[i])));

                temp = temp + f[i];
                c[i] = temp;
            }
        }

        //eta : update root
        vector<double>eta(n-1);
        I1.clear();I1.resize(0);I1.shrink_to_fit();
        I2.clear();I2.resize(0);I2.shrink_to_fit();

        for(int i = 0; i < a.size(); ++i){
            if(a[i]>0)
                I1.push_back(i);
            else
                I2.push_back(i);
        }

        for(int i = 0; i < I1.size(); ++i){
            eta[I1[i]] = 2*b[I1[i]] / (a[I1[i]] + std::sqrt(std::abs((a[I1[i]] * a[I1[i]]) - 4 * b[I1[i]]*c[I1[i]])));                         
        }

        for(int i = 0; i < I2.size(); ++i ){
            eta[I2[i]] = (a[I2[i]] - std::sqrt(std::abs((a[I2[i]] * a[I2[i]]) - 4*b[I2[i]]*c[I2[i]]))) / (2*c[I2[i]]);            
        }

        vector<double> Ic;
        for(int i = 0 ; i < c.size(); i++){
            if(std::abs(c[i]) == 0)
                Ic.push_back(i);
        }

        if(Ic.size() != 0){
            vector<double> Ia;
            for(int i = 0; i < Ic.size(); i++){
                if(std::abs(a[Ic[i]]) == 0)
                    Ia.push_back(i);
            }

            if(Ia.size() != 0){
                I_B.clear();I_B.resize(0);I_B.shrink_to_fit();
                for(int i = 0; i < Ia.size(); ++i){
                    double temp = Ic[Ia[i]];
                    I_B.push_back(temp);
                }

                for(int i = 0; i< I_B.size(); i++){
                    double temp;
                    if(f0[(int)I_B[i]] >= 0){
                        temp = tempvSqr[(int)I_B[i]] + ((d[(int)I_B[i]+1] - x[(int)I_B[i]])*(d[(int)I_B[i]+1] - x[(int)I_B[i]])) * ((df[(int)I_B[i]] - tempvSqr[(int)I_B[i]]) / ((d[(int)I_B[i]] - x[(int)I_B[i]])*(d[(int)I_B[i]] - x[(int)I_B[i]])));                         
                    }
                    else{
                        temp = tempvSqr[(int)I_B[i]+1] + std::pow((d[(int)I_B[i]] - x[(int)I_B[i]]), 2.0) * (df[(int)I_B[i]] - tempvSqr[(int)I_B[i]+1] / std::pow((d[(int)I_B[i]+1] - x[(int)I_B[i]]),2.0));
                    }
                    a[(int)I_B[i]] = temp;
                }             
            }

            for(int i = 0; i < Ic.size(); ++i){
                eta[(int)Ic[i]] = b[int(Ic[i])] / a[(int)Ic[i]];
            }
        }

        // f*eta should be negative, otherwise run a newton step
        I_B.clear(); I_B.resize(0); I_B.shrink_to_fit();

        for(int i = 0; i < (n-1); i++){
            if(f[i]*eta[i] >= 0)
                I_B.push_back(i);
        }

        if(I_B.size() != 0){
            for(int i = 0 ; i < I_B.size(); i++){
                eta[(int)I_B[i]] = -f[(int)I_B[i]] / df[(int)I_B[i]];
            }
        }

        //bound check with bisection safeguard
        vector<double> tmp;
        for(int i = 0; i < tau.size(); ++i){
            double temp = tau[i] + eta[i];
            tmp.push_back(temp);
        }

        //  I = find( tmp<xlb | tmp>xub | tmp==0 | tmp.*xub<0 | tmp.*xlb<0 | isnan(tmp) );
        I_B.clear(); I_B.resize(0); I_B.shrink_to_fit();
        for(int i = 0; i < tmp.size(); ++i){
            if(tmp[i] < xlb[i] | tmp[i] > xub[i] | tmp[i] == 0 | tmp[i]*xub[i] < 0 | tmp[i]*xlb[i] < 0 | std::isnan(tmp[i])){
                I_B.push_back(i);
            }
        }

        for(int i = 0; i < I_B.size(); ++i){
            double temp;
            if(f[(int)I_B[i]] >= 0)
                temp = xlb[(int)I_B[i]];
            else
                temp = xub[(int)I_B[i]];

            temp = temp - tau[(int)I_B[i]];
            eta[(int)I_B[i]] = temp / 2;

        }                   


        //update root 
        for(int i = 0; i < tau.size(); ++i)
            tau[i] = tau[i] + eta[i];
        
        for(int i = 0; i < org.size(); i++){
            x[i] = tau[i] + d[(int)org[i]];
        }

        iter_ct = iter_ct + 1;
    //end of while loop
    }
   //TODO: handle the memory allocation/deallocation for f, psi, phi etc.

    // check residual after MAX_ITER iterations
    
    if(FMM_ITER){
        if(n >=N){
            //trifmm1dlocal shift


        //end of inner if block    
        }
        else{
            int kRows = org.size();
            int kCols  = d.size();             
            double *K = new double[kRows*kCols];

            for(int row = 0; row < kRows; ++row){
                for(int col=0; col < kCols; ++col){
                    K[col+row*(kCols)] = d[row] - d[col];
                }
            }

            double* temp_tau = &tau[0];

            //double* temp_org;
            //memcpy(temp_org,temp_o,sizeof(int)*org.size());
            bsxfun('P',&K,{kRows,kCols},temp_tau,{tau.size(),1});

            for(int i = 0; i < (kRows*kCols); ++i)
                K[i] = 1/K[i];

            double *K1 = new double[kRows*kCols];
            double *K2 = new double[kRows*kCols];

            memcpy(K1,K,sizeof(double)*kRows*kCols);
            memcpy(K2,K,sizeof(double)*kRows*kCols);

            //K1 = tril(K), K2 = triu(K,1)
            for(int row = 0; row < kRows; row++){
                for(int col = 0; col < kCols; col++){
                    if(row < col)
                        K1[col + row*(kCols)] = 0;
                    else
                        K2[col + row*(kCols)] = 0;
                }
            }

            psi = new double[kRows];
            phi = new double[kRows];

            memset(psi,0,sizeof(double)*kRows);
            memset(phi,0,sizeof(double)*kRows);

            double *temp = new double[v.size()];
            for(int i = 0; i < v.size(); ++i)
                temp[i] = v[i]*v[i];

            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,kRows,1,kCols,1,K1,kCols,temp,1,0,psi,1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,kRows,1,kCols,1,K2,kCols,temp,1,0,phi,1);
            
            f = new double[kRows];

            for(int i = 0; i < kRows; i++)
                f[i] = rho-psi[i]-phi[i];

            pd_size = kRows;
         //end of inner else block   
        }
        
        vector<double> res;
        residual.clear(); residual.resize(0); residual.shrink_to_fit();
        J2.clear(); J2.resize(0); J2.shrink_to_fit();
        for(int i = 0; i < pd_size; i++){

            residual.push_back(std::abs(f[i]));
            if(std::abs(f[i]) > C*n*eps*std::max(res[i],alpha))
                J2.push_back(i);

        }   
        //add print statements for continuing iterations 
    //end of outer if block    
    }

    else{
         J2.clear(); J2.resize(0); J2.shrink_to_fit();
         for(int i = 0; i < (n-1); ++i)
            J2.push_back(i);

    //end of outer else block    
    }

    // Find root x(J2)
    

//not completed still working on it, don't use this routine
}
