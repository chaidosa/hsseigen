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

double dot(vector<double>& a, vector<double>& b, int n){
    double result = 0;
    for(int i = 0; i < n; ++i)
        result +=a[i]*b[i];

    return result;
}

double PSI_1(vector<double>& v, double x, int iter, vector<double>& delta){
    vector<double>temp_v;
    vector<double>temp_delta;
    for(int i = 0; i < iter; i++){
        temp_v.push_back((v[i]*v[i]));
        double temp = delta[i] - x;
        temp_delta.push_back(1/temp);
    }
    return dot(temp_v,temp_delta, iter);
}

double DPSI_1(vector<double>& v, double x, int iter, vector<double>& delta){
    vector<double>temp_v;
    vector<double>temp_delta;
    for(int i = 0; i < iter; i++){
        temp_v.push_back((v[i]*v[i]));

        double temp = delta[i] - x;
        temp = temp*temp;
        temp_delta.push_back(1/temp);
    }
    return dot(temp_v,temp_delta, iter);
}

double PSI_2(vector<double>&v, double x, vector<double>&delta){
    vector<double>temp_v;
    vector<double>temp_d;
    for(int i = 0; i < v.size()-1; i++){
        temp_v.push_back(v[i]*v[i]);

        double temp = delta[i]-x;
        temp_d.push_back(1/temp);
    }
    return dot(temp_v, temp_d, temp_v.size());
}

double DPSI_2(vector<double>&v, double x, vector<double>&delta){
    vector<double>temp_v;
    vector<double>temp_d;
    for(int i = 0; i < v.size()-1; i++){
        temp_v.push_back(v[i]*v[i]);

        double temp = delta[i]-x;
        temp = temp*temp;
        temp_d.push_back(1/temp);
    }
    return dot(temp_v, temp_d, temp_v.size());
}

double PHI_2(vector<double>&v, double x, vector<double>& delta){
    double temp = v[v.size()-1];
    temp = temp *temp;
    temp = temp / (delta[delta.size()-1] - x);

    return temp;
}

double DPHI_2(vector<double>&v, double x, vector<double>& delta){
    double temp = v[v.size()-1];
    temp = temp *temp;
    double temp2 =  (delta[delta.size()-1] - x);
    temp2 = temp2 * temp2;

    return temp/temp2;
}


double PHI_1(vector<double>& v, double x, int iter, int n, vector<double>& delta){
    
    vector<double> temp_v;
    vector<double> temp_delta;

    for(int i =n-1; i > i; --i){        
        temp_v.push_back((v[i]*v[i]));
        
        double temp = delta[i] - x;
        temp_delta.push_back(1/temp);        
    }

    return dot(temp_v,temp_delta, temp_v.size());
}

double DPHI_1(vector<double>& v, double x, int iter, int n, vector<double>& delta){

    vector<double> temp_v;
    vector<double> temp_delta;

    for(int i =n-1; i > i; --i){        
        temp_v.push_back((v[i]*v[i]));
        
        double temp = delta[i] - x;
        temp = temp*temp;
        temp_delta.push_back(1/temp);        
    }

    return dot(temp_v,temp_delta, temp_v.size());

}


void printArray(double **Arr, int row, int col){
    ofstream txtOut;
    txtOut.open("output.txt", std::ofstream::out | std::ofstream::app);
    double *A = *Arr;
    for(int i = 0; i < row; i++){
        for(int j=0; j < col; j++){
            txtOut<<A[j+i*col]<<"\t";
        }
        txtOut<<"\n";
    }
    txtOut<<"\n";
    txtOut.close();
}


void rootfinder(vector<double>& d,vector<double>& v){
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
        std::vector<double>x;
        std::vector<double>tau;
        std::vector<double>org;
        double temp = (v[0]*v[0]) + d[0];
        x.push_back(temp);
        temp = temp - d[0];
        tau.push_back(temp);
        org.push_back(0);
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
    int stop_criteria = 1;
    int display_warning =1;
   // if(record){
        std::vector<double> iter(n);
        std::vector<double> residual(n);
        std::vector<double> erretms(n);
  //  }
    double v_norm = vec_norm(v);
    double rho    = 1 / (v_norm*v_norm);

    for(int i = 0; i < v.size(); ++i){
        v[i] = v[i] / v_norm;
    }
    std::vector<double> d0 = diff_vec(d);

    std::vector<double> x;
    for(int i = 0; i < n-1; ++i){
        double temp = (d[i] + d[i+1]) / 2;
        x.push_back(temp);
    }

    std::vector<double> org;
    for(int i=0; i < n-1; ++i){
        org.push_back(i);
    }

   double *f0;
   int f0_size = 0;
   //sqare of vector v
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
       
        int kRows = org.size();
        int kCols  = d.size(); 
        
        double *K = new double[kRows*kCols];

        for(int row = 0; row < kRows; ++row){
            for(int col=0; col < kCols; ++col){
                K[col+row*(kCols)] = d[(int)org[row]] - d[col];
            }
        }

        //printArray(&K, kRows, kCols);
        
        //it used to hold d0/2
        double *tempD = new double[d0.size()];

        for(int i = 0; i < d0.size(); ++i){
            tempD[i] = d0[i] / 2;
        }

        bsxfun('P', &K, {kRows,kCols}, tempD, {d0.size(), 1});
        
        //  K = 1 ./ K;
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
        delete[] tempD;
        K = NULL;
    }

    // h = 2 * diff(v.^2) ./ d0;
    std::vector<double>tempvSqr(v2_arr, v2_arr+v.size());
    std::vector<double> h = diff_vec(tempvSqr);
    for(int i = 0; i < h.size(); ++i){
        h[i] = 2*h[i] / d0[i];
    }

    //g = f0 - h;
    std::vector<double>g;
    for(int i = 0; i < f0_size; ++i){
        double temp = (f0[i]-h[i]);
        g.push_back(temp);
    }

    // org = [1:n-1]' + (f0<0);
    for(int i = 0; i < org.size(); ++i){
        if(f0[i] < 0){
            org[i] = org[i] + 1;
        }
    }




    // **Upper bound and lower bound**
    std::vector<double> d1(d.begin(), d.end()-1);
    std::vector<double> d2(d.begin()+1, d.end());
    std::vector<double> v1(v.begin(), v.end()-1);
    std::vector<double> v2(v.begin()+1, v.end());
    std::vector<double> xub(n-1,0);
    std::vector<double> xlb(n-1,0);
    std::vector<int> I1;
    std::vector<int> I2;

    //find (f0 >= 0) and find(f0 < 0)
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

    std::vector<double>xlb0(xlb.begin(), xlb.end());
    std::vector<double>xub0(xub.begin(), xub.end());




    // **inital guess**      
    std::vector<double>a; //a = (2*(f0>=0) - 1) .* d0 .* g + v(1:n-1).^2 + v(2:n).^2;
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
    std::vector<double>b;
    for(int i = 0; i < n-1; ++i){
        double temp_first_part = 0;
        double temp_second_part = 0; 

        if(f0[i] >= 0)
            temp_first_part = tempvSqr[i];          
        else            
            temp_second_part = tempvSqr[i+1];        

        temp_first_part = (temp_first_part - temp_second_part)*d0[i];
        b.push_back(temp_first_part);
    }

    std::vector<double>tau(n-1, 0);    
    //clearnig previous allocations
    I1.resize(0); I1.shrink_to_fit(); I2.resize(0); I2.shrink_to_fit();

    for(int i = 0; i < a.size(); ++i){

        if(a[i]>0)
            I1.push_back(i);    
        else
            I2.push_back(i);
        
    }

    //tau(I1) = 2*b(I1) ./ (a(I1) + sqrt(abs(a(I1).^2 - 4*b(I1).*g(I1))));    
    for(int i = 0; i < I1.size(); ++i){
        tau[I1[i]] = (2*b[I1[i]]) / (a[I1[i]] + std::sqrt(std::abs( (a[I1[i]]*a[I1[i]]) - (4*b[I1[i]]) * (g[I1[i]]))));
    }

    //tau(I2) = (a(I2) - sqrt(abs(a(I2).^2 - 4*b(I2).*g(I2)))) ./ (2*g(I2));
    for(int i = 0; i < I2.size(); ++i){
        tau[I2[i]] = (a[I2[i]] - std::sqrt(std::abs( (a[I2[i]]* a[I2[i]]) - (4 * b[I2[i]]*g[I2[i]])))) / (2* g[I2[i]]);
    }




    // **bound check**
    //this vector is same as I in matlab version
    std::vector<int>I_B;
    for(int i = 0; i < tau.size(); ++i){

        if((tau[i]<=xlb[i]) | (tau[i]>=xub[i]))
            I_B.push_back(i);

    }

    if(I_B.size() != 0)
        for(int i = 0; i < I_B.size(); ++i)
            tau[I_B[i]] = (xub[I_B[i]] + xlb[I_B[i]])/2;      
    

    for(int i = 0; i < tau.size(); ++i)
        x[i] = tau[i] + d[(int)org[i]];




    
    /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%  iterate n-1 roots simultaneously  %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    int iter_ct = 0; //iteration count
    std::vector<double> swtch(n-1, 1);
    std::vector<double> pref;
    std::vector<double> f_prev;
    std::vector<double> J2;
    double *f;
    double  *psi, *phi, *dpsi, *dphi, *df;     
        int pd_size = 0;

    while (iter_ct < FMM_ITER && Flag){    
        if(iter_ct > 0)
            pref = f_prev;

    
        // **bound check with bisection safeguard**
        I_B.clear(); I_B.resize(0); I_B.shrink_to_fit();
        for(int i = 0; i < tau.size(); i++){
            if(tau[i] == 0 | (tau[i]*xub[i]) < 0 | (tau[i]*xlb[i])<0 | isnan(tau[i]) | tau[i] < xlb0[i] | tau[i] > xub0[i])
                I_B.push_back(i);
        }

        if(~I_B.empty()){
            for(int i = 0; i < I_B.size(); ++i){
                tau[I_B[i]] = (xub[I_B[i]] + xlb[I_B[i]])/2;
                x[I_B[i]] = tau[I_B[i]] + d[(int)org[I_B[i]]];
            }
        }




        // **secular function evaluation**     
        if(n >= N){
            //fmm
        }

        else{

            int kRows = org.size();
            int kCols  = d.size();             
            double *K = new double[kRows*kCols];

            for(int row = 0; row < kRows; ++row){
                for(int col=0; col < kCols; ++col){
                    K[col+row*(kCols)] = d[(int)org[row]] - d[col];
                }
            }
            
            bsxfun('P', &K, {kRows,kCols}, &tau[0], {tau.size(), 1});

            for(int i = 0; i < (kRows*kCols); ++i)
                K[i] = 1/K[i];

            double *K1 = new double[kRows*kCols];
            double *K2 = new double[kRows*kCols];

            memcpy(K1, K, sizeof(double)*kRows*kCols);
            memcpy(K2, K, sizeof(double)*kRows*kCols);

            //K1 = tril(K), K2 = triu(K,1)
            for(int row = 0; row < kRows; row++){
                for(int col = 0; col < kCols; col++){
                    if(row < col)
                        K1[col + row*(kCols)] = 0;
                    else
                        K2[col + row*(kCols)] = 0;
                }
            }
           // printArray(&K1, kRows, kCols);            
           // printArray(&K2, kRows, kCols);

            psi = new double[kRows];
            phi = new double[kRows];

            memset(psi, 0, sizeof(double)*kRows);
            memset(phi, 0, sizeof(double)*kRows);            

            //psi = K1*v^2; phi = K2*v^2;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1, K1, kCols, v2_arr, 1, 0, psi, 1);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1, K2, kCols, v2_arr, 1, 0, phi, 1);
            
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

            memset(dpsi, 0, sizeof(double)*kRows);
            memset(dphi, 0, sizeof(double)*kRows);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1, dK1, kCols, v2_arr, 1, 0, dpsi, 1);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1, dK2, kCols, v2_arr, 1, 0, dphi, 1);

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



        // **adaptive FMM iterations**        
        std::vector<double> res(pd_size);

        for(int i = 0 ; i < pd_size; ++i)
            res[i] = rho + std::abs(psi[i]) + std::abs(phi[i]);

        double nJ2 = 0;
        if(iter_ct > 0)
            nJ2 = J2.size();        
        
        for(int i = 0; i < pd_size; ++i)       
            if(std::abs(f[i]) > (C*n*eps*std::max(res[i], alpha)))
                J2.push_back(i);         
        
                
        if(iter_ct > 0 && (J2.size() <= 0.01*n || std::abs(nJ2 - J2.size()) < 0.01 *n))
            Flag = Flag - 1;

        if(iter_ct == 5 || iter_ct < 5 && Flag == 0)
            percent = J2.size() / (n-1);
        


        // **update  upper and lower bounds**
        I_B.clear(); I_B.resize(0); I_B.shrink_to_fit();        
        std::vector<int> II;
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



        // **swtch = 1 : fixed weight method**
        // **swtch = 0 : middle way method**
        std::vector<int> Iswtch;
        if(iter_ct > 0){
            for(int i = 0; i < pd_size; ++i){
                if((pref[i]*f[i] > 0) & (std::abs(pref[i]) > abs(f[i]) / 10))
                    Iswtch.push_back(i);
            }
            for(int i = 0; i < Iswtch.size(); ++i){
                swtch[Iswtch[i]] = -swtch[Iswtch[i]] + 1;
            }
        }



        // **quadratic equation coefficients**
        a.clear(); a.resize(0); a.shrink_to_fit();
        b.clear(); b.resize(0); b.shrink_to_fit();
        // a = ((d1 - x) + (d2 - x)) .* f - (d1 - x) .* (d2 - x) .* df;
        // b = (d1 - x) .* (d2 - x) .* f;
        for(int i = 0; i < pd_size; ++i){
            double temp_a = ( (d[i]-x[i]) + (d2[i]-x[i]) ) * f[i] - (d1[i]-x[i]) * (d2[i]-x[i]) * df[i];
            a.push_back(temp_a);

            double temp_b = (d1[i]-x[i]) * (d2[i]-x[i]) * f[i];
            b.push_back(temp_b);
        }

        std::vector<double>c(n-1, 0);
        if(iter_ct > 0){
            std::vector<double> c0;
            std::vector<double> c1;
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



        // **eta : update root**
        std::vector<double>eta(n-1);
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

        std::vector<double> Ic; //handle corner case ?
        for(int i = 0 ; i < c.size(); i++){
            if(std::abs(c[i]) == 0)
                Ic.push_back(i);
        }

        if(Ic.size() != 0){

            std::vector<double> Ia;
            for(int i = 0; i < Ic.size(); i++){
                if(std::abs(a[Ic[i]]) == 0)
                    Ia.push_back(i);
            }

            if(Ia.size() != 0){
                I_B.clear(); I_B.resize(0); I_B.shrink_to_fit();
                for(int i = 0; i < Ia.size(); ++i){
                    double temp = Ic[(int)Ia[i]];
                    I_B.push_back(temp);
                }

                for(int i = 0; i< I_B.size(); i++){
                    double temp;
                    if(f0[(int)I_B[i]] >= 0){
                        temp = tempvSqr[I_B[i]] + ((d[I_B[i]+1] - x[I_B[i]])*(d[I_B[i]+1] - x[I_B[i]])) * ((df[I_B[i]] - tempvSqr[I_B[i]]) / ((d[(int)I_B[i]] - x[(int)I_B[i]])*(d[(int)I_B[i]] - x[(int)I_B[i]])));                         
                    }
                    else{
                        temp = tempvSqr[I_B[i]+1] + std::pow((d[(int)I_B[i]] - x[(int)I_B[i]]), 2.0) * (df[(int)I_B[i]] - tempvSqr[(int)I_B[i]+1] / std::pow((d[(int)I_B[i]+1] - x[(int)I_B[i]]),2.0));
                    }
                    a[I_B[i]] = temp;
                }             
            }

            for(int i = 0; i < Ic.size(); ++i){
                eta[(int)Ic[i]] = b[int(Ic[i])] / a[(int)Ic[i]];
            }
        }



        // **f*eta should be negative, otherwise run a newton step**
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



        //**bound check with bisection safeguard**
        std::vector<double> tmp;
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



        // **update root** 
        for(int i = 0; i < tau.size(); ++i)
            tau[i] = tau[i] + eta[i];
        
        for(int i = 0; i < org.size(); i++){
            x[i] = tau[i] + d[(int)org[i]];
        }

        iter_ct = iter_ct + 1;
        f_prev.insert(f_prev.begin(),f,f+pd_size);
    //end of while loop
        delete[] psi, phi, dpsi, dphi, f;
    }
   //TODO: handle the memory allocation/deallocation for f, psi, phi etc.




    // **check residual after MAX_ITER iterations**    
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
                    K[col+row*(kCols)] = d[(int)org[row]] - d[col];
                }
            }          

            bsxfun('P',&K,{kRows,kCols},&tau[0],{tau.size(),1});

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

            delete[] K, K1, K2;
            delete[] psi, phi;
         //end of inner else block   
        }
        
        std::vector<double> res;
        for(int i = 0; i < pd_size; ++i)
            res.push_back(rho + std::abs(psi[i]) + std::abs(phi[i]));

        residual.clear(); residual.resize(0); residual.shrink_to_fit();
        J2.clear(); J2.resize(0); J2.shrink_to_fit();
        for(int i = 0; i < pd_size; i++){

            residual.push_back(std::abs(f[i]));
            if(std::abs(f[i]) > C*n*eps*std::max(res[i],alpha))
                J2.push_back(i);

        }
        if(n > N0);   
            //add print statements for continuing iterations 
    //end of outer if block    
    }

    else{
        if(n > N0);
            //add print statements for continuing iterations
        J2.clear(); J2.resize(0); J2.shrink_to_fit();

        for(int i = 0; i < (n-1); ++i)
            J2.push_back(i);

    //end of outer else block  
    
    }




/*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% find roots x(J2) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
    std::vector<double> QJ(n, 0);
    
    for(int k = 0; k < J2.size(); k++){
        int i = J2[k];
        //double y0 = (d[i] + d[i+1]) / 2;        
        //double w  = rho + PSI_1(v, y0, i, d) + PHI_1(v, y0, i, n, d);
        double y0, w;
        double org0;

        vector<double>delta;        

        double ub,lb = 0;       
        double dw, erretm;
        double x0;
    /*    
        if(w >= 0){
            org0 = i; //origin
            for(int iter = 0; iter < d.size(); iter++) //shift
                delta.push_back(d[iter] - d[i]);
            
            ub = delta[i+1] / 2;
            lb = 0;

            if(lb < tau[i] && tau[i] <= ub)
                x0 = tau[i];
            else
                x0 = (ub + lb) / 2;              

        }
        else if (w < 0){
            org0 = i + 1; //origin
            for(int iter = 0; iter < d.size(); iter++) //shift
                delta.push_back(d[iter] - d[i+1]);
            
            ub = 0;
            lb = delta[i] / 2;

            if(lb < tau[i] && tau[i] < ub)
                x0 = tau[i];
            else
                x0 = (ub + lb) / 2;
        }
    */    
        //** shift origin according to the newly created w0 **       

        double x0 = d[i] / 2;

        for(int j = 0; j < d.size(); j++)
            delta.push_back(d[j] - d[i]);
        double w0 = rho + PSI_1(v, x0, i, delta);

        if(w0 >= 0){
            org0 = i;
            org[i] = org0;
            ub = delta[i+1] / 2;
            lb = 0;
            if((lb < tau[i]) && (tau[i] <= ub))
                x0 = tau[i];
            else
                x0 = (ub + lb) / 2;

        }
        else if(w0 < 0){
            org0 = i + 1;
            org[i] = org0;
            delta.clear(); delta.resize(0); delta.shrink_to_fit();
            for(int iter = 0; iter < d.size(); iter++) 
                delta.push_back(d[iter] - d[i+1]);
            ub = 0;
            lb = delta[i] / 2;
            if((lb < tau[i]) && (tau[i] < ub))
                x0 = tau[i];
            else
                x0 = (ub + lb) / 2;
        }



        //** shift origin according to the newly created f0(i) **
        /*
        %     org0 = org(i);      
        %     delta = d - d(org0);
        %     x0 = tau(i);
        %     ub = xub(i);
        %     lb = xlb(i);
        */


        //**initial lower and upper bound**
        double ub0 = ub;
        double lb0 = lb;



        //**rational inerpolation**
        double psi0 = PSI_1(v, x0, i, delta);
        double phi0 = PHI_1(v, x0, i, v.size(), delta);
        double dpsi0 = DPSI_1(v, x0, i, delta);
        double dphi0 = DPHI_1(v, x0, i, v.size(), delta);

        w = rho + (psi0 + phi0);
        dw = dpsi0 + dphi0;
        if(stop_criteria == 0){

        }
        else  if(stop_criteria == 1)
            erretm = C * n * (rho + std::abs(psi0) + std::abs(phi0));

        iter_ct = 0;
        double swtch_2 = 1;

        while(std::abs(w) > erretm * eps){
            double prew;

            // **maximal iterations**
            if (iter_ct >= MAX_ITER){
                if(display_warning){
                    printf("root does not converge");
                    break;
                }
            }


            //**update upper and lower bounds**
            if (w >= 0)
                ub = std::min(ub, x0);
            else
                lb = std::max(lb, x0);


            
            //swtch_2 = 1 : fixed weight method
            //swtch_2 = 0 : middle way method
            if (iter_ct > 0)
                if((prew * w > 0) && (std::abs(prew) > std::abs(w) / 10))
                    swtch_2 = -swtch_2 + 1;


            
            //Quadratic equation coefficients
            double A = (delta[i] - x0 + delta[i+1] -x0) * w - (delta[i] - x0)*(delta[i+1] - x0)*dw;
            double B = (delta[i] - x0) * (delta[i+1] - x0)*w;
            double c;
            double eta;
            if(swtch_2 == 0)
                c = - (delta[i+1] - x0) * DPSI_1(v,x0,i,delta) - (delta[i+1] - x0) * DPHI_1(v, x0, i, v.size(), delta) + w ;
            else if (swtch_2 == 1){
                if(f0[i] >= 0)
                    c = -(delta[i+1] - x0) * dw - (v[i]*v[i]) *(delta[i] -delta[i+1])/((delta[i] - x0)*(delta[i] - x0)) + w;
                else if(f0[i] < 0)
                    c = -(delta[i] - x0) * dw - (v[i+1]*v[i+1]) * (delta[i+1] - delta[i]) / ((delta[i+1] - x0)*(delta[i+1] - x0)) + w;    
            }

            //eta :: **root update**
            if (abs(c) == 0){  // handle corner case : c==0
                if(abs(A) == 0){
                    if(f0[i] >= 0)
                        A = (v[i] * v[i]) + ((delta[i+1] - x0)*(delta[i+1] - x0)) * (dw - (v[i]*v[i]) / ((delta[i] - x0) * (delta[i] - x0)));
                    else if (f0[i] < 0)
                        A = (v[i+1]*v[i+1]) + ((delta[i] - x0)*(delta[i] - x0)) * (dw - (v[i+1]*v[i+1]) /((delta[i+1] - x0) *(delta[i+1] - x0)));
                }
                eta = B / A;   
            }
            else{
                if(A <= 0)
                    eta = (A - std::sqrt(std::abs((A*A) - 4*B*C))) / (2*c);
                else
                    eta = 2*B / (A + std::sqrt(std::abs(A*A - 4*B*c)));
            }

            // ** f*eta should be negative, otherwise run a newton step **
            if(w*eta >= 0)
                eta = -w / dw;
            
            // **check if updated root lies in the [lb, ub]**
            if(((x0 + eta) < lb) || ((x0 + eta) > ub) || ((x0 + eta) <= lb0) || ((x0 + eta) >= ub0)){
                if(w >= 0)
                    eta = (lb - x0) / 2;
                else
                    eta = (ub - x0) / 2;
            }

            //**update root**
            x0 = x0 + eta;



            //**for next iteration**
            iter_ct = iter_ct + 1;
            prew  = w;
            psi0  = PSI_1(v, x0, i, delta);
            phi0  = PHI_1(v, x0, i, v.size(), delta);
            dpsi0 = DPSI_1(v, x0, i, delta);
            dphi0 = DPHI_1(v, x0, i, v.size(), delta);
                w = rho + psi0 + phi0;
               dw = dpsi0 + dphi0;

            if(stop_criteria == 0){
                //later   
            } 
            else if(stop_criteria == 1)
                erretm = C * n* (rho + std::abs(psi0) + std::abs(phi0));               
                   
        }

        if(record){
            residual[i] = std::abs(w);
            erretms[i]  = std::abs(erretm);
            iter[i]     = iter_ct;
        }


        // i^th root and eigenvector
        if(org0 == i){
            tau[i] = x0;
            x[i]   = x0 + d[i];
        }
        else if(org0 == i+1){
            tau[i] = x0;
            x[i]   = x0 + d[i+1];
        }        
    }



/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% find n^th root %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
   
   //delta = d - d(n);
    vector<double> delta;
    for(int iter = 0; iter < d.size(); iter++){
        double temp = d[iter] - d[d.size()-1];
        delta.push_back(temp);
    }


    //initial guess
    double x0 = 1 / (2 * rho);
    double lb0  = 0;
    double w = rho + PSI_2(v, x0, delta) + PHI_2(v, x0, delta);
    double c = w - (v[n-2]*v[n-2]) / (delta[n-2] - x0) - (v[n-1]*v[n-1]) / (delta[n-1] - x0);       
    double A, B, lb, ub,dw,erretm,eta;

    if(w <= 0){
        double temp = (v[n-2]*v[n-2]) / (1 / rho - delta[n-2]) + (v[n-1]*v[n-1])*rho;
        if(c <= temp)
            x0 = 1 / rho;
        else{
            A = c * delta[n-2] + (v[n-2]*v[n-2]) + (v[n-1]*v[n-1])/(delta[n-1] - x0);
            B = -(v[n-1]*v[n-1])*delta[n-2];

            if(A < 0)
                x0 = 2*B / (std::sqrt(std::abs(A*A + 4*B*c)) - A); 
            else
                x0 = (A + std::sqrt(std::abs(A*A + 4*B*c))) / (2*c);
        }
        lb = 1 / (2 * rho);
        ub = 1 / rho;
    }
    else{
        A = c*delta[n-2] + (v[n-2]*v[n-2]) + (v[n-1]*v[n-1]);        
        B = - (v[n-1]*v[n-1]) * delta[n-2];
        if(A < 0)
            x0 = 2*B / (std::sqrt(std::abs(A*A + 4*B*c))-A);
        else
            x0 = (A + std::sqrt(std::abs(A*A + 4*B*c))) / (2*c);

        lb = 0;
        ub = 1 / (2*rho);
    }



    //**rational interpolation**
    double psi0 = PSI_2(v, x0, delta);
    double phi0 = PHI_2(v, x0, delta);
    double dpsi0 = DPSI_2(v,x0,delta);
    double dphi0 = DPHI_2(v, x0, delta);
    w = rho + psi0 + phi0;
    dw = dpsi0 + dphi0;

    if(stop_criteria == 0){
        //later
    }
    else if(stop_criteria == 1)
        erretm = C * n* (rho + std::abs(psi0) + std::abs(phi0));
    
    iter_ct = 0;
    while (std::abs(w) > erretm*eps)
    {
        //**maximum iterations**
        if(iter_ct >= MAX_ITER){
            if(display_warning)
                printf("root doesn't converge");
            break;
        }


        //**update upper and lower bounds**
        if(w >= 0)
            ub = std::min(ub, x0);
        else
            lb = std::max(lb, x0);
        


        //calculate new root 
        A = (delta[n-2] - x0 + delta[n-1] - x0) * w - (delta[n-2] - x0) * (delta[n-1] - x0)* dw;
        B = (delta[n-2] - x0) * (delta[n-1]-x0) *w;
        c = -(delta[n-2] - x0) * DPSI_2(v,x0,delta) - (delta[n-1] - x0) * DPHI_2(v,x0,delta) + w;
        c = std::abs(c);
        if(c == 0)
            eta = ub - x0;
        else{
            if(A >= 0)
                eta = (A + std::sqrt(std::abs(A*A - 4*B*C))) / (2*c);
            else
                eta = 2*B / (A - std::sqrt(std::abs(A*A - 4*B*C)));
        }



        // ** f*eta should be negative, otherwise run a Newton step**
        if(w*eta >= 0)
            eta = -w / dw;


        
        // **check if updated root lies in the [lb, ub]**
        if((x0 + eta) < lb || ((x0 + eta) > ub) || ((x0 + eta) <= lb0)){
             if(w >= 0)
                eta = (lb-x0) / 2;
            else
                eta = (ub-x0) / 2;
        }
        

        //**update root**
        x0 = x0 + eta;


        // **for next iteration**
        psi0 = PSI_2(v, x0, delta);
        phi0 = PHI_2(v, x0, delta);
        dpsi0 = DPSI_2(v, x0, delta);
        dphi0 = DPHI_2(v, x0, delta);
        w = rho + psi0 + phi0;
        dw = dpsi0 + dphi0;

        if(stop_criteria == 0){
            //later
        }
        else if(stop_criteria == 1)
            erretm = C * n* (rho + std::abs(psi0) + std::abs(phi0));

        iter_ct = iter_ct + 1;
    }


    if(record){
        residual[n-1] = std::abs(w);
        erretms[n-1]  = std::abs(erretm);
        iter[n-1] = iter_ct;
    }
    

    //**n^th root and eigenvector**
    tau.push_back(x0);
    org.push_back(n-1);
    x.push_back(x0 + d[n-1]);

}