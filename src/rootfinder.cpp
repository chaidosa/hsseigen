#include <bits/stdc++.h>
#include <iostream>
#include <string>
#include <cassert>
#include "eigenmatrix.h"
#include "bsxfun.h"
#include "rootfinder.h"
#include "fmm1d_local_shift_2.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}

using namespace std;

double dot(vector<double>& a, vector<double>& b, int n)
{
    double result = 0;
    for(unsigned int i = 0; i <= (unsigned)n; ++i)
        result +=a[i]*b[i];

    return result;
}


double PSI_1(vector<double>& v, double x, int iter, vector<double>&delta)
{
    vector<double>temp_v;
    vector<double>temp_delta;
    for(unsigned int i = 0; i <=(unsigned)iter; i++)
    {
        temp_v.push_back((v[i]*v[i]));
        double temp = delta[i] - x;
        temp_delta.push_back(1/temp);
    }
    return dot(temp_v,temp_delta, iter);
}


double DPSI_1(vector<double>& v, double x, int iter, vector<double>& delta)
{
    vector<double>temp_v;
    vector<double>temp_delta;
    for(unsigned int i = 0; i <=(unsigned)iter; i++)
    {
        temp_v.push_back((v[i]*v[i]));

        double temp = delta[i] - x;
        temp = temp*temp;
        temp_delta.push_back(1/temp);
    }
    return dot(temp_v,temp_delta, iter);
}


double PSI_2(vector<double>&v, double x, vector<double>&delta)
{
    vector<double>temp_v;
    vector<double>temp_d;
    for(unsigned int i = 0; i <(unsigned)v.size()-1; i++)
    {
        temp_v.push_back(v[i]*v[i]);

        double temp = delta[i]-x;
        temp_d.push_back(1/temp);
    }
    return dot(temp_v, temp_d, temp_v.size()-1);
}


double DPSI_2(vector<double>&v, double x, vector<double>&delta)
{
    vector<double>temp_v;
    vector<double>temp_d;
    for(unsigned int i = 0; i <(unsigned)v.size()-1; i++)
    {
        temp_v.push_back(v[i]*v[i]);

        double temp = delta[i]-x;
        temp = temp*temp;
        temp_d.push_back(1/temp);
    }
    return dot(temp_v, temp_d, temp_v.size()-1);
}


double PHI_2(vector<double>&v, double x, vector<double>& delta)
{
    double temp = v[v.size()-1];
    temp = temp *temp;
    temp = temp / (delta[delta.size()-1] - x);

    return temp;
}


double DPHI_2(vector<double>&v, double x, vector<double>& delta)
{
    double temp = v[v.size()-1];
    temp = temp *temp;
    double temp2 =  (delta[delta.size()-1] - x);
    temp2 = temp2 * temp2;

    return temp/temp2;
}


double PHI_1(vector<double>& v, double x, int iter, int n, vector<double>& delta)
{
    
    vector<double> temp_v;
    vector<double> temp_delta;

    for(unsigned int i =n-1; i >(unsigned)iter; --i)
    {        
        temp_v.push_back((v[i]*v[i]));
        
        double temp = delta[i] - x;
        temp_delta.push_back(1/temp);        
    }

    return dot(temp_v,temp_delta, temp_v.size()-1);
}


double DPHI_1(vector<double>& v, double x, int iter, int n, vector<double>& delta)
{

    vector<double> temp_v;
    vector<double> temp_delta;

    for(unsigned int i =n-1; i >(unsigned)iter; --i)
    {        
        temp_v.push_back((v[i]*v[i]));
        
        double temp = delta[i] - x;
        temp = temp*temp;
        temp_delta.push_back(1/temp);        
    }

    return dot(temp_v,temp_delta, temp_v.size()-1);

}
 

// This for testing the output of matrix
void printArray(double **Arr, int row, int col)
{
    ofstream txtOut;
    txtOut.open("output.txt", std::ofstream::out | std::ofstream::app);
    double *A = *Arr;
    for(unsigned int i = 0; i <(unsigned)row; i++)
    {
        for(int j=0; j < col; j++)
        {
            txtOut<<A[j+i*col]<<"\t";
        }
        txtOut<<"\n";
    }
    txtOut<<"\n";
    txtOut.close();
}


Root *rootfinder(vector<double>& d,vector<double>& v, double N)
{
    /*
    %%% Input:
    %%% d, v: as in secular equation
    %%% N: size threshold to use fmm

    %%% Ouput:
    %%% x: roots of secular equation 
    %%% tau: gaps
    %%% org: shifts for each root
    */


    Root * results;

    //double N = 17000;    
    int dSize = d.size();
    int n = v.size();
    if(dSize != n)
    {
        cout << "not computable\n";
        assert(false);
    }
    //edge case
    if(n == 1)
    {
	    results = new Root(1);
	    double temp = (v[0]*v[0]) + d[0];
	    results->x=new double(temp);
            results->org = new int(0);
            temp = temp - d[0];
            results->tau.push_back(temp);        
            results->percent = 1;
            return results;
    }

    results = new Root(n-1);
    
    double C = 64;
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
   
    std::vector<double> iter(n, 0);
    std::vector<double> residual(n, 0);
    std::vector<double> erretms(n, 0);
 
    double v_norm = vec_norm(v);
    double rho    = 1 / (v_norm*v_norm);

    results->x=new double[n];
    double* x=results->x;

    int *org = new int[n];
    int org_size = n-1;

    double *v2_arr;   
    v2_arr = new double[n];  

    double *tempD = new double[dSize - 1];  

    std::vector<double>d0(n-1);
    for(int i = 0; i < n-1; ++i)
    {
        v[i] = v[i] / v_norm;
        
        d0[i] = d[i+1] - d[i];
        
        x[i] = (d[i] + d[i+1]) / 2;
        
        org[i] = i;

        v2_arr[i] = v[i]*v[i];  
    
        tempD[i] = d0[i] / 2;    
    
    }
    v[n-1] = v[n-1] / v_norm;
    v2_arr[n-1] = v[n-1]*v[n-1];

    //std::vector<double> d0 = diff_vec(d);

    //results->x=new double[n];
    //double* x=results->x;
    /*
    for(unsigned int i = 0; i <(unsigned)n-1; ++i)
    {
        double temp = (d[i] + d[i+1]) / 2;
        x[i]=temp;
    }
    */
    //std::vector<int> org;
   /* int *org = new int[n];
    int org_size = n-1;
    for(unsigned int i=0; i <(unsigned) n-1; ++i)
        org[i] = i;
    */

   int kRows = org_size;
   int kCols  = dSize; 
   double *f0;
   f0 = new double[kRows]; 
   int f0_size = kRows;
   
   
   //sqares of vector v so that v^2 isn't computed again
   /*
   double *v2_arr;   
   v2_arr = new double[n];
   for(unsigned int i = 0; i < n; ++i)
        v2_arr[i] = v[i]*v[i];
    */

   //it used to hold d0/2
   //double *tempD = new double[dSize - 1];

   //for(unsigned int i = 0; i <(unsigned)dSize - 1; ++i)
   // tempD[i] = d0[i] / 2;
   
   if(n >= N)
   {
	   
        double* z = trifmm1d_local_shift(r, x, d.data(), v2_arr, tempD, org, 1, org_size, dSize, 1);
        //f0 = rho - fl - fu;
        for(unsigned int i = 0; i <(unsigned)kRows; ++i)
            f0[i] = rho - z[i] - z[kRows+i];
   }
   else
   {
       //  K = d(org) - d.';   
        double *K = new double[kRows*kCols];

        for(unsigned int row = 0; row <(unsigned)kRows; ++row)
            for(unsigned int col=0; col <(unsigned)kCols; ++col)
                K[col+row*(kCols)] = d[org[row]] - d[col];          
        
        //  K = bsxfun(@plus, K, d0/2);       
        bsxfun('P', &K, make_pair(kRows,kCols), tempD, make_pair(dSize - 1, 1));        
        
        //  K = 1 ./ K;
        for(unsigned int row_col = 0; row_col <(unsigned)(kRows*kCols); ++row_col)
            K[row_col] = 1 / K[row_col];
        
        memset(f0, 0, sizeof(double)*kRows);
        // f0 = rho - K * v.^2;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1.0, K, kCols, v2_arr, 1, 0.0, f0, 1);

        for(unsigned int i = 0; i <(unsigned)kRows; ++i)
            f0[i] = rho - f0[i];
        
        delete[] K;
        K = NULL;
    }
    
    delete[] tempD;

    // h = 2 * diff(v.^2) ./ d0;

    //can be made efficient will do in next round
    std::vector<double>tempvSqr(v2_arr, v2_arr+ n);
    //std::vector<double> h = diff_vec(tempvSqr);
    std::vector<double> h(n-1);
    std::vector<double>g(f0_size);


    for(int i = 0; i < n-1; ++i)
    {
        h[i] = 2*(v2_arr[i+1] - v2_arr[i])/ d0[i];

        g[i] = f0[i] - h[i];

        if(f0[i] < 0)
            org[i] = org[i] + 1;
    }
   
    
    
    // **Upper bound and lower bound**
    //std::vector<double> d1(d.begin(), d.end()-1); //only referenced no need to copy
        double *d1 = d.data();
        double *d2 = d.data() + 1;
    //std::vector<double> d2(d.begin()+1, d.end()); //only referenced
    
    //std::vector<double> v1(v.begin(), v.end()-1); //only referenced
        double *v1 = v.data();
        double *v2 = v.data() + 1;
    //std::vector<double> v2(v.begin()+1, v.end()); // only referenced
    std::vector<double> xub(n-1,0);
    std::vector<double> xlb(n-1,0);
    std::vector<int> I1;
    I1.reserve(n-1); // reserving more will find an efficient way later
    std::vector<int> I2;
    I2.reserve(n-1);
    //find (f0 >= 0) and find(f0 < 0)
    for(unsigned int i = 0; i <(unsigned)f0_size; i++)
    {
        if(f0[i]>=0)
            I1.push_back(i);
        else
            I2.push_back(i);
    }

    for(unsigned int i = 0; i <(unsigned)I1.size(); ++i)
        xub[I1[i]] = (d2[I1[i]] - d1[I1[i]]) / 2;  
    

    for(unsigned int i = 0; i <(unsigned)I2.size(); ++i)
        xlb[I2[i]] = (d1[I2[i]] - d2[I2[i]]) / 2;  
    

    //std::vector<double>xlb0(xlb); //only referenced
    //std::vector<double>xub0(xub);
    double *xlb0 = xlb.data();
    double *xub0 = xub.data();

    // **inital guess**      
    
    std::vector<double>a(n-1); //a = (2*(f0>=0) - 1) .* d0 .* g + v(1:n-1).^2 + v(2:n).^2;
    std::vector<double>b(n -1); //b = ((f0>=0) .* v(1:n-1).^2 - (f0<0) .* v(2:n).^2) .* d0;

    for(unsigned int i = 0; i <(unsigned)n-1; i++)
    {
        double temp_first_part = 0;
        double temp_second_part = 0;
    
        temp_first_part = f0[i]>=0?1:-1;

        temp_first_part *= d0[i]*g[i];
        temp_second_part = tempvSqr[i]+tempvSqr[i+1];
        temp_first_part += temp_second_part;        
        a[i] = temp_first_part;

        temp_first_part = 0;
        temp_second_part = 0; 

        if(f0[i] >= 0)
            temp_first_part = tempvSqr[i];          
        else            
            temp_second_part = tempvSqr[i+1];        

        temp_first_part = (temp_first_part - temp_second_part)*d0[i];     
        b[i] = temp_first_part;

    }    
    

    std::vector<double>tau(org_size, 0);    
    //clearnig previous allocations
    I1.clear(); I2.clear(); 

    for(unsigned int i = 0; i <(unsigned)a.size(); i++)
    {

        if(a[i]>0)
            I1.push_back(i);    
        else
            I2.push_back(i);
        
    }

    //tau(I1) = 2*b(I1) ./ (a(I1) + sqrt(abs(a(I1).^2 - 4*b(I1).*g(I1))));    
    for(unsigned int i = 0; i <(unsigned)I1.size(); i++)
    {
        tau[I1[i]] = (2*b[I1[i]]) / (a[I1[i]] + std::sqrt(std::abs( (a[I1[i]]*a[I1[i]]) - (4*b[I1[i]]) * (g[I1[i]]))));
    }

    //tau(I2) = (a(I2) - sqrt(abs(a(I2).^2 - 4*b(I2).*g(I2)))) ./ (2*g(I2));
    for(unsigned int i = 0; i <(unsigned)I2.size(); i++)
    {
        tau[I2[i]] = (a[I2[i]] - std::sqrt(std::abs( (a[I2[i]]* a[I2[i]]) - (4 * b[I2[i]]*g[I2[i]]) ))) / (2* g[I2[i]]);
    }


    // **bound check**
    //this vector is same as I in matlab version
    std::vector<int>I_B;
    I_B.reserve(org_size);
    for(unsigned int i = 0; i <(unsigned)tau.size(); i++)
        if((tau[i]<=xlb[i]) | (tau[i]>=xub[i]))
            I_B.push_back(i);

    if(!I_B.empty())
        for(unsigned int i = 0; i <(unsigned)I_B.size(); ++i)
            tau[I_B[i]] = (xub[I_B[i]] + xlb[I_B[i]])/2;      
    

    for(unsigned int i = 0; i <(unsigned)tau.size(); ++i)
        x[i] = tau[i] + d[org[i]];



    
    /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%  iterate n-1 roots simultaneously  %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    int iter_ct = 0; //iteration count
    std::vector<double> swtch(n-1, 1);
    std::vector<double> pref;
    std::vector<double> f_prev;
    std::vector<double> J2;
    double *f = new double[org_size];

    double  *psi=NULL, *phi=NULL, *dpsi=NULL, *dphi=NULL, *df=NULL;
    double *K=NULL, *K1=NULL, *K2=NULL;

    int pd_size = 0;
    if(n < N){
        df = new double[kRows];            
        dpsi = new double[kRows];
        psi = new double[kRows];
        dphi = new double[kRows];        
        phi = new double[kRows];
        K = new double[kRows*kCols];
        K1 = new double[kRows*kCols];
        K2 = new double[kRows*kCols];  
    }

    while ((iter_ct < FMM_ITER) && Flag)
    {    
        if(iter_ct > 0)
        {
            pref.clear();
            copy(f_prev.begin(), f_prev.end(), back_inserter(pref));            
        }
    
        // **bound check with bisection safeguard**
        I_B.clear();
        for(unsigned int i = 0; i <(unsigned)tau.size(); i++)
        {
            if((tau[i] == 0) | ((tau[i]*xub[i]) < 0) | ((tau[i]*xlb[i]) < 0) | (std::isnan(tau[i])) | (tau[i] < xlb0[i]) | (tau[i] > xub0[i]))
                I_B.push_back(i);
        }

        if(!I_B.empty())
        {
            for(unsigned int i = 0; i <(unsigned)I_B.size(); ++i)
            {
                tau[I_B[i]] = (xub[I_B[i]] + xlb[I_B[i]]) / 2;
                x[I_B[i]] = tau[I_B[i]] + d[org[I_B[i]]];
            }
        }



        kRows = org_size;
        kCols = dSize;             
       // f = new double[kRows];
    /*
        dpsi = new double[kRows];
        dphi = new double[kRows];
        df = new double[kRows];

        memset(dpsi, 0, sizeof(double)*kRows);
        memset(dphi, 0, sizeof(double)*kRows);
        psi = new double[kRows];
        phi = new double[kRows];

    */

        // **secular function evaluation**     
        if(n >= N)
        {
            //fmm
	        double* z = trifmm1d_local_shift(r, x, d.data(), v2_arr, tau.data(), org, 1, org_size, dSize, 1);
            for(unsigned int i = 0; i <(unsigned)kRows; ++i)
              f[i] = rho - z[i] - z[kRows+i];
            psi = z; phi = z+kRows; 

            double *zd = trifmm1d_local_shift(r, x, d.data(), v2_arr, tau.data(), org, 1, org_size, dSize, 2);
            dpsi = zd; dphi = zd+kRows;
            
            for(unsigned int i = 0; i <(unsigned)kRows; ++i)
              df[i] = zd[i] + zd[kRows+i];
        }

        else
        {
                   
            memset(dpsi, 0, sizeof(double)*kRows);
            memset(dphi, 0, sizeof(double)*kRows);           

            
            memset(psi, 0, sizeof(double)*kRows);
            memset(phi, 0, sizeof(double)*kRows);  
            
            memset(K, 0, sizeof(double)*kRows*kCols);

            for(unsigned int row = 0; row <(unsigned)kRows; ++row)
            {
                for(unsigned int col = 0; col <(unsigned)kCols; ++col)
                {
                    K[col+row*(kCols)] = d[org[row]] - d[col];
                }
            }
            
            bsxfun('P', &K, {kRows,kCols}, tau.data(), {tau.size(), 1});            

            for(unsigned int i = 0; i <(unsigned)(kRows*kCols); ++i)
                K[i] = 1/K[i];            

            memcpy(K1, K, sizeof(double)*kRows*kCols);
            memcpy(K2, K, sizeof(double)*kRows*kCols);

            //K1 = tril(K), K2 = triu(K,1)
            for(unsigned int row = 0; row <(unsigned)kRows; row++)
            {
                for(unsigned int col = 0; col <(unsigned)kCols; col++)
                {
                    if(row < col)
                        K1[col + row*(kCols)] = 0;
                    else
                        K2[col + row*(kCols)] = 0;
                }
            }
           // printArray(&K1, kRows, kCols);            
           // printArray(&K2, kRows, kCols);


                    

            //psi = K1*v^2; phi = K2*v^2;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1, K1, kCols, v2_arr, 1, 0, psi, 1);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1, K2, kCols, v2_arr, 1, 0, phi, 1);
            

            for(unsigned int i = 0; i <(unsigned)kRows; i++)
                f[i] = rho-psi[i]-phi[i];
            
            double *dK1 = K1;
            double *dK2 = K2;

            for(unsigned int row = 0; row <(unsigned)kRows; row++)
            {
                for(unsigned int col = 0; col <(unsigned)kCols; col++)
                {

                    dK2[col + row*kCols] = dK2[col + row*kCols]*dK2[col + row*kCols];        
                    dK1[col + row*kCols] = dK1[col + row*kCols]*dK1[col + row*kCols];

                }
            }


            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1, dK1, kCols, v2_arr, 1, 0, dpsi, 1);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kRows, 1, kCols, 1, dK2, kCols, v2_arr, 1, 0, dphi, 1);

            
            for(unsigned int i = 0; i <(unsigned)kRows; ++i)
            {
                df[i] = dpsi[i] + dphi[i];
            }

            pd_size = kRows;           
                    
            
        // end of else statement
        }



        // **adaptive FMM iterations**        
        std::vector<double> res(pd_size);

        for(unsigned int i = 0 ; i <(unsigned)pd_size; ++i)
            res[i] = rho + std::abs(psi[i]) + std::abs(phi[i]);

        double nJ2 = 0;
        if(iter_ct > 0)
            nJ2 = J2.size();        
        
        J2.clear();
        for(unsigned int i = 0; i <(unsigned)pd_size; ++i)       
            if(std::abs(f[i]) > (C*n*eps*(std::max(res[i], alpha))))
                J2.push_back(i);         
        
                
        if(iter_ct > 0 && (J2.size() <= 0.01*n || std::abs(nJ2 - J2.size()) < 0.01 *n))
            Flag = Flag - 1;

        if(iter_ct == 5 || (iter_ct < 5 && Flag == 0))
            percent = J2.size() / (n-1);
        


        // **update  upper and lower bounds**
        I_B.clear();     
        std::vector<int> II;
        I1.reserve(pd_size);
        for(unsigned int i = 0; i <(unsigned)pd_size; i++)
        {            
            if(f[i] < 0){
                //I_B.push_back(i);
                xlb[i] = std::max(tau[i], xlb[i]);
            }
            else{
                //II.push_back(i);
                xub[i] = std::min(tau[i], xub[i]);
            }
        }

        // **swtch = 1 : fixed weight method**
        // **swtch = 0 : middle way method**
        std::vector<int> Iswtch;
        Iswtch.reserve(pd_size);

        if(iter_ct > 0)
        {
            for(unsigned int i = 0; i <(unsigned)pd_size; ++i)
            {
                if((pref[i]*f[i] > 0) & (std::abs(pref[i]) > abs(f[i]) / 10))
                    Iswtch.push_back(i);
            }

            for(unsigned int i = 0; i <(unsigned)Iswtch.size(); ++i)
            {
                swtch[Iswtch[i]] = -swtch[Iswtch[i]] + 1;
            }

        }



        // **quadratic equation coefficients**
        a.clear();
        b.clear(); 
        // a = ((d1 - x) + (d2 - x)) .* f - (d1 - x) .* (d2 - x) .* df;
        // b = (d1 - x) .* (d2 - x) .* f;
        for(unsigned int i = 0; i <(unsigned)pd_size; ++i)
        {
            double temp_var = ( (d1[i]-x[i]) + (d2[i]-x[i]) ) * f[i] - (d1[i]-x[i]) * (d2[i]-x[i]) * df[i];
            a.push_back(temp_var);

            temp_var = (d1[i]-x[i]) * (d2[i]-x[i]) * f[i];
            b.push_back(temp_var);
        }

        std::vector<double>c(n-1, 0);
        if(iter_ct > 0)
        {
            std::vector<double> c0(pd_size);
            std::vector<double> c1(pd_size);
            for(unsigned int i = 0; i < (unsigned)pd_size; ++i)
            {

                double temp = -(d1[i]-x[i]) * dpsi[i] - (d2[i]-x[i]) * dphi[i] + f[i];
                c0[i] = temp;

                temp = 0;
                if(f0[i] >= 0)
                    temp = -1*( (d2[i] - x[i]) * df[i] + (v1[i]*v1[i]) * (d1[i]- d2[i]) / std::pow((d1[i]-x[i]), 2));
                else
                    temp = -1*((d1[i] - x[i]) * df[i] + (v2[i]*v2[i])*(d2[i]-d1[i]) / std::pow((d2[i]-x[i]), 2));

                temp = temp + f[i];
                c1[i] =temp;

            }
            
            for(unsigned int i = 0; i <(unsigned)swtch.size(); ++i)
            {
                if(swtch[i] == 0)
                    c[i] = c0[i];
                else
                    c[i] = c1[i];

            }
        }

        else
        {
            for(unsigned int i = 0; i <(unsigned)pd_size; ++i)
            {
                double temp;
                if(f0[i]>=0)
                    temp = -1*( (d2[i] - x[i]) * df[i] + (v1[i]*v1[i]) * (d1[i]- d2[i]) / ( (d1[i]-x[i]) * (d1[i]-x[i]) ));
                else
                    temp = -1*( (d1[i] - x[i]) * df[i] + (v2[i]*v2[i])*(d2[i]-d1[i]) / ((d2[i]-x[i]) * (d2[i]-x[i])));

                temp = temp + f[i];
                c[i] = temp;
            }
        }



        // **eta : update root**
        std::vector<double>eta(n-1, 0);
        I1.clear();
        I2.clear();
        for(unsigned int i = 0; i < a.size(); ++i)
        {
            if(a[i]>0){
                I1.push_back(i);
                eta[i] = (2*b[i]) / ( a[i] + std::sqrt(std::abs( (a[i] * a[i]) - 4 * b[i]*c[i]) ) );   
            }
            else{
                I2.push_back(i);
                eta[i] = (a[i] - std::sqrt(std::abs((a[i] * a[i]) - 4*b[i]*c[i]))) / (2*c[i]);
            }
        }    
        
        std::vector<int> Ic; //handle corner case ?
        Ic.reserve(c.size());
        for(unsigned int i = 0 ; i <(unsigned)c.size(); i++)
        {
            if(std::abs(c[i]) == 0)
                Ic.push_back(i);
        }

        if(Ic.size() != 0)
        {

            std::vector<int> Ia;
            Ia.reserve(Ic.size());
            for(unsigned int i = 0; i <(unsigned)Ic.size(); i++)
            {
                if(std::abs(a[Ic[i]]) == 0)
                    Ia.push_back(i);
            }

            if(Ia.size() != 0)
            {
                I_B.clear();
                for(unsigned int i = 0; i <(unsigned)Ia.size(); ++i)
                {
                    double temp = Ic[Ia[i]];
                    I_B.push_back(temp);
                }

                for(unsigned int i = 0; i<(unsigned)I_B.size(); i++)
                {
                    double temp;
                    if(f0[I_B[i]] >= 0)
                    {
                        temp = tempvSqr[I_B[i]] + std::pow((d[I_B[i]+1] - x[I_B[i]]), 2) * (df[I_B[i]] - tempvSqr[I_B[i]] / (std::pow((d[I_B[i]] - x[I_B[i]]), 2)));                         
                    }
                    else
                    {
                        temp = tempvSqr[I_B[i]+1] + std::pow((d[I_B[i]] - x[I_B[i]]), 2.0) * (df[I_B[i]] - tempvSqr[I_B[i]+1] / std::pow((d[I_B[i]+1] - x[I_B[i]]),2.0));
                    }
                    a[I_B[i]] = temp;
                }             
            }

            for(unsigned int i = 0; i <(unsigned)Ic.size(); ++i)
            {
                eta[Ic[i]] = b[Ic[i]] / a[Ic[i]];
            }
        }



        // **f*eta should be negative, otherwise run a newton step**        

        for(unsigned int i = 0; i <(unsigned)(n-1); i++)
        {
            if(f[i]*eta[i] >= 0){               
                eta[i] = -f[i] / df[i];
            }
        }       


        //**bound check with bisection safeguard**
        std::vector<double> tmp;
        for(unsigned int i = 0; i <(unsigned)tau.size(); ++i){
            double temp = tau[i] + eta[i];
            tmp.push_back(temp);
        }

        //  I = find( tmp<xlb | tmp>xub | tmp==0 | tmp.*xub<0 | tmp.*xlb<0 | isnan(tmp) );
        
        for(unsigned int i = 0; i <(unsigned)tmp.size(); ++i)
        {
            if((tmp[i] < xlb[i]) | (tmp[i] > xub[i]) | (tmp[i] == 0) | (tmp[i]*xub[i] < 0) | (tmp[i]*xlb[i] < 0) | (std::isnan(tmp[i])))
            {                
                double temp = (f[i] >= 0) ? xlb[i] : xub[i];        
                temp = temp - tau[i];
                eta[i] = temp / 2;
            }
        }  

        // **update root** 
        for(unsigned int i = 0; i <(unsigned)tau.size(); i++){
            tau[i] = tau[i] + eta[i];
            x[i] = tau[i] + d[org[i]];
        }
        
   
        iter_ct = iter_ct + 1;
        f_prev.insert(f_prev.begin(),f,f+pd_size);  
        //end of while loop
    }
   
    delete[] psi;	
	delete[] dpsi;
    delete[] f;
    delete[] df;
	psi=NULL;
	phi=NULL;
	dpsi=NULL;
	dphi=NULL;
	f=NULL;
	df=NULL; 
    if(n < N){
        delete [] K;
        delete [] K1;
        delete [] K2;
        delete [] phi;
        delete[] dphi;
    }


    // **check residual after MAX_ITER iterations**
    //psi = new double[kRows];
    //phi = new double[kRows];

    //memset(psi,0,sizeof(double)*kRows);
    //memset(phi,0,sizeof(double)*kRows);
    f = new double[kRows];

    if(FMM_ITER)
    {
        if(n >=N){
            //trifmm1dlocal shift
            double * z = trifmm1d_local_shift(r, x, d.data(), v2_arr, tau.data(), org, 1, org_size, dSize, 1);
            for(unsigned int i = 0; i <(unsigned)kRows; ++i)
              f[i] = rho - z[i] - z[kRows+i];       
           
            psi = z; phi = z+kRows;
        }
        else{
            int kRows = org_size;
            int kCols  = dSize;             
            double *K = new double[kRows*kCols];

            for(unsigned int row = 0; row <(unsigned)kRows; ++row){
                for(unsigned int col=0; col <(unsigned)kCols; ++col){
                    K[col+row*(kCols)] = d[org[row]] - d[col];
                }
            }    
            
            bsxfun('P',&K,{kRows,kCols},&tau[0],{tau.size(),1});
        
            for(unsigned int i = 0; i <(unsigned)(kRows*kCols); ++i)
                K[i] = 1/K[i];

            double *K1 = new double[kRows*kCols];
            double *K2 = new double[kRows*kCols];

            memcpy(K1,K,sizeof(double)*kRows*kCols);
            memcpy(K2,K,sizeof(double)*kRows*kCols);

            //K1 = tril(K), K2 = triu(K,1)
            for(unsigned int row = 0; row <(unsigned)kRows; row++){
                for(unsigned int col = 0; col <(unsigned)kCols; col++){
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

            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,kRows,1,kCols,1,K1,kCols,v2_arr,1,0,psi,1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,kRows,1,kCols,1,K2,kCols,v2_arr,1,0,phi,1);            
            

            for(unsigned int i = 0; i <(unsigned)kRows; i++)
                f[i] = rho-psi[i]-phi[i];

            pd_size = kRows;

            delete[] K;
	        delete[] K1;
	        delete[] K2;
         //end of inner else block   
        }
        
        residual.clear(); 
        J2.clear();
        std::vector<double> res(pd_size);
        for(unsigned int i = 0; i <(unsigned)pd_size; ++i){
            res[i] = (rho + std::abs(psi[i]) + std::abs(phi[i]));
            residual.push_back(std::abs(f[i]));
            if(std::abs(f[i]) > C*n*eps*std::max(res[i],alpha))
                J2.push_back(i);
        }
     
        if(n > N0)
            printf("Roofinder: continue to iterate");
          
            
        delete[] psi;
	    delete[] phi;
    //end of outer if block    
    }

    else{
        if(n > N0)
            printf("rootfinder: continue to iterate");
            //add print statements for continuing iterations
        J2.clear();
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
    
    for(unsigned int k = 0; k < J2.size(); k++)
    {
        int current_i = J2[k];
      
        double w;
        double org0;

        std::vector<double>delta;        
        delta.reserve(dSize);
        double ub,lb = 0;       
        double dw, erretm;
        double x0;

        x0 = d0[current_i] / 2;

        for(unsigned int j = 0; j <(unsigned)dSize; j++)
            delta.push_back(d[j] - d[current_i]);

        double w0 = rho + PSI_1(v, x0, current_i, delta) + PHI_1(v, x0, current_i, v.size(), delta);

        if(w0 >= 0){
            org0 = current_i;
            org[current_i] = org0;
            ub = delta[current_i+1] / 2;
            lb = 0;
            if((lb < tau[current_i]) && (tau[current_i] <= ub))
                x0 = tau[current_i];
            else
                x0 = (ub + lb) / 2;

        }
        else if(w0 < 0){
            org0 = current_i + 1;
            org[current_i] = org0;
            delta.clear(); 
            for(unsigned int iter = 0; iter <(unsigned)dSize; iter++) 
                delta.push_back(d[iter] - d[current_i+1]);
            ub = 0;
            lb = delta[current_i] / 2;
            if((lb < tau[current_i]) && (tau[current_i] < ub))
                x0 = tau[current_i];
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
        double psi0 = PSI_1(v, x0, current_i, delta);
        double phi0 = PHI_1(v, x0, current_i, n, delta);
        double dpsi0 = DPSI_1(v, x0, current_i, delta);
        double dphi0 = DPHI_1(v, x0, current_i, n, delta);

        w = rho + (psi0 + phi0);
        dw = dpsi0 + dphi0;

        if(stop_criteria == 0);
             
        else  if(stop_criteria == 1)
            erretm = C * n * (rho + std::abs(psi0) + std::abs(phi0));

        iter_ct = 0;
        double swtch_2 = 1;
        
        double prew;
        while(std::abs(w) > erretm * eps){
            

            // **maximal iterations**
            if (iter_ct >= MAX_ITER){
                if(display_warning){
                    printf("root does not converge");
                    break;
                }
            }


            //**update upper and lower bounds**
            if (w >= 0)
                ub = ub > x0 ? x0 : ub;                
            else
                lb = lb > x0 ? lb : x0;               


            
            //swtch_2 = 1 : fixed weight method
            //swtch_2 = 0 : middle way method
            if (iter_ct > 0)
                if((prew * w > 0) && (std::abs(prew) > std::abs(w) / 10))
                    swtch_2 = -swtch_2 + 1;


            
            //Quadratic equation coefficients
            double A = (delta[current_i] - x0 + delta[current_i+1] -x0) * w - (delta[current_i] - x0)*(delta[current_i+1] - x0) * dw;
            double B = (delta[current_i] - x0) * (delta[current_i+1] - x0) * w;
            double c;
            double eta;
            if(swtch_2 == 0)
                c = - (delta[current_i] - x0) * DPSI_1(v,x0,current_i,delta) - (delta[current_i+1] - x0) * DPHI_1(v, x0, current_i, n, delta) + w ;
            else if (swtch_2 == 1)
            {
                if(org0 == current_i)
                    c = -(delta[current_i+1] - x0) * dw - (tempvSqr[current_i]) * (delta[current_i] - delta[current_i+1]) / std::pow((delta[current_i] - x0), 2) + w;
                else if(org0 == current_i + 1)
                    c = -(delta[current_i] - x0) * dw - (tempvSqr[current_i + 1]) * (delta[current_i+1] - delta[current_i]) / std::pow((delta[current_i+1] - x0), 2) + w;    
            }

            //eta :: **root update**
            if (abs(c) == 0){  // handle corner case : c==0
                if(abs(A) == 0){
                    if(org0 == current_i)
                        A = (tempvSqr[current_i]) + std::pow((delta[current_i+1] - x0), 2) * (dw - tempvSqr[current_i] / std::pow((delta[current_i] - x0), 2));
                    else if (org0 == current_i + 1)
                        A = (tempvSqr[current_i+1]) + std::pow((delta[current_i] - x0), 2) * (dw - (tempvSqr[current_i + 1]) / std::pow((delta[current_i+1] - x0), 2));
                }
                eta = B / A;   
            }
            else{

                if(A <= 0)
                    eta = (A - std::sqrt(std::abs((A*A) - 4*B*c))) / (2*c);
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
            psi0  = PSI_1(v, x0, current_i, delta);
            phi0  = PHI_1(v, x0, current_i, n, delta);
            dpsi0 = DPSI_1(v, x0, current_i, delta);
            dphi0 = DPHI_1(v, x0, current_i, n, delta);
                
            w = rho + psi0 + phi0;
            dw = dpsi0 + dphi0;

            if(stop_criteria == 0){
                //later   
            } 
            else if(stop_criteria == 1)
                erretm = C * n* (rho + std::abs(psi0) + std::abs(phi0));               
                   
        }

        if(record){
            residual[current_i] = std::abs(w);
            erretms[current_i]  = std::abs(erretm);
            iter[current_i]     = iter_ct;
        }


        // i^th root and eigenvector
        if(org0 == current_i){
            tau[current_i] = x0;
            x[current_i]   = x0 + d[current_i];
        }
        else if(org0 == current_i+1){
            tau[current_i] = x0;
            x[current_i]   = x0 + d[current_i+1];
        }        
    }



/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% find n^th root %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
   
   //delta = d - d(n);
    vector<double> delta;
    for(unsigned int iter = 0; iter < dSize; iter++){
        double temp = d[iter] - d[dSize - 1];
        delta.push_back(temp);
    }


    //initial guess
    double x0 = 1 / (2 * rho);
    double lb0  = 0;
    double w = rho + PSI_2(v, x0, delta) + PHI_2(v, x0, delta);
    double c = w - (tempvSqr[n-2]) / (delta[n-2] - x0) - (tempvSqr[n-1]) / (delta[n-1] - x0);       
    double A, B, lb, ub,dw,erretm,eta;

    if(w <= 0){
        double temp = (tempvSqr[n-2]) / (1 / rho - delta[n-2]) + (tempvSqr[n-1])*rho;
        if(c <= temp)
            x0 = 1 / rho;
        
        else{
            A = c * delta[n-2] + (tempvSqr[n-2]) + (tempvSqr[n-1]);
            B = -(tempvSqr[n-1])*delta[n-2];

            if(A < 0)
                x0 = 2*B / (std::sqrt(std::abs(A*A + 4*B*c)) - A); 
            else
                x0 = (A + std::sqrt(std::abs(A*A + 4*B*c))) / (2*c);
        }

        lb = 1 / (2 * rho);
        ub = 1 / rho;
    }
    else
    {
        A = c*delta[n-2] + (tempvSqr[n-2]) + (tempvSqr[n-1]);        
        B = - (tempvSqr[n-1]) * delta[n-2];
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
                eta = 2*B / (A - std::sqrt(std::abs(A*A - 4*B*c)));
        }



        // ** f*eta should be negative, otherwise run a Newton step**
        if(w*eta >= 0)
            eta = -w / dw;


        
        // **check if updated root lies in the [lb, ub]**
        if((x0 + eta) < lb || ((x0 + eta) > ub) || ((x0 + eta) <= lb0))
        {
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


   // if(record){
        residual[n-1] = std::abs(w);
        erretms[n-1]  = std::abs(erretm);
        iter[n-1] = iter_ct;
    //}
    
    delete [] v2_arr;
    delete [] f0;
    delete [] f;
    //**n^th root and eigenvector**
    tau.push_back(x0);
    org[n-1] = n-1;
    org_size++;;
    x[n-1]=x0 + d[n-1];
    results->tau = tau;
    results->org = org;
    results->org_size = org_size;
    results->x   = x;
    results->percent = percent;
    return results;
}
