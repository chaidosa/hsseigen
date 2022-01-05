#include "secular.h"
#include<bits/stdc++.h>
#include <assert.h>
#include <vector>
#include <set>
#include<algorithm>
#include<iostream>
#include "eigenmatrix.h"

using namespace std;

//norm of a vector using eucledian method
double vec_norm(vector<double> v){
    double result = 0;
    for(int i = 0; i < v.size(); ++i){
        result +=v[i]*v[i];
    }

return sqrt(result);    
}

//Function similar to matlab: calculates diff of next and current element and form a vector
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
    vector<double>x;
    vector<double>tau;
    vector<double>org;

    int n = v.size();
    if(d.size()!=v.size()){
        cout << "not computable";
        assert(false);
    }
    //edge case
    if(n == 1){
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

    /*
    if n >= N
    [fl, fu, nflops1] = trifmm1d_local_shift(r, x, d, v.^2, d0/2, org, 1);
    f0 = rho - fl - fu;

else
    K = d(org) - d.';
    K = bsxfun(@plus, K, d0/2);
    K = 1 ./ K;
    f0 = rho - K * v.^2;
    nflops = nflops + flops('prod', K, 'n', v, 'n');
end
h = 2 * diff(v.^2) ./ d0;
g = f0 - h;
org = [1:n-1]' + (f0<0);

    */
   if(n >= N){
        //[fl, fu, nflops1] = trifmm1d_local_shift(r, x, d, v.^2, d0/2, org, 1);
        //f0 = rho - fl - fu;
   }
   else{    

   }

}

nonleaf* secular(double *d, int dSize, double *v, int vSize,double N=1024){
/*    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% roots of secular equation w/ deflation %%%%%
%%%%% 1 + \sum_j (v_j^2 / (lam1_j - lam)) = 0%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input:
% d,v: as in secular equation
% tol: tolerance for deflation (will add later)
% N:   size threshold to use fmm (will add later)
%
%%% Ouput:
% Lam: roots of secular equation
% Q:   structured eigen matrix
*/

    if(dSize != vSize){
        std::cout<<"diagonal and vector size should be same\n";
        assert(false);
    }
    //to store the results of secular 
    nonleaf *result = new nonleaf();

    double tol = 1.0e-10;
    int n      = vSize;
//step 1: deflate small vi
    std::vector<int> Tempt;
    std::vector<int> zero_2_n;
    for(int i = 0; i < vSize; ++i){
        if(std::abs(v[i]) < tol){
            Tempt.push_back(i);
        }
        //creating another array having indices 1:n
        zero_2_n.push_back(i);
    }

    std::vector<int>diff;
    std::set_difference(zero_2_n.begin(),zero_2_n.end(),Tempt.begin(),Tempt.end(),std::inserter(diff, diff.begin()));

    std::vector<double> Lam1;
    for(int i = 0 ; i < Tempt.size() ; ++i){
        Lam1.push_back(d[Tempt[i]]);
    }
    std::vector<double> Lam2;
    std::vector<double> Lam3;

    std::vector<double>d2;
    for(int i = 0; i < diff.size(); ++i){
        d2.push_back(d[diff[i]]);
    }

    std::vector<double>v2;
    for(int i = 0; i < diff.size(); ++i){
        v2.push_back(v[diff[i]]);
    }

    int n1 = Tempt.size();
    for(int i = 0; i < diff.size(); ++i){
        Tempt.push_back(diff[i]);
    }

// step 2: deflate  close eigenvalue
    if(n1 < n){
       //conj(v2) ./ abs(v2); 
       std::vector<double>v2c;
       for(int i = 0; i < v2.size(); ++i){
           int temp = v2[i] / abs(v2[i]);
           v2c.push_back(temp);
           v2[i] = abs(v2[i]);
       }
       
       //sort array by keeping track of previous index
       
       //vector to store element       
        vector<pair<double,int>> vP;
        for(int i = 0; i < d2.size(); ++i){
            vP.push_back(make_pair(d[i],i));
        }
        //sorting pair vector
        sort(vP.begin(),vP.end());

        vector<double> tempV2;
        tempV2 = v2;
        for(int i = 0 ;i < v2.size(); ++i){
            int index = vP[i].second;
            v2[i] = tempV2[index];
        }

        int p=1;
        vector<double> J,G;
        for(int j = 1; j < d2.size(); ++j){
            double s   = v2[p];
            double c   = v2[j];
            double tau = sqrt((c*c)+(s*s));
            double t   = d2[j] - d2[p];
            c          = c / tau;
            s          = s / tau;

            // close eigenvalues, deflate d2(p) 
            if(abs(t*c*s) < tol){
                v2[p] = 0;
                v2[j] = tau;
                t     = ((c*c)*(d2[p])) + ((s*s)*(d2[j]));
                d2[j] = (s*s)*(d2[p]) + (c*c)*(d2[j]);
                d2[p] = t;
                // Givens rotation
                G.push_back(p);
                G.push_back(j);
                G.push_back(c);
                G.push_back(s);
                J.push_back(p);
                p = j;
            }
            else{
                p = j;
            }
        }
   
        std::vector<int>Jc;
        vector<int> n_sub_n1;
        for(int i = 0; i < (n-n1); ++i){
            n_sub_n1.push_back(i);
        }
        std::set_difference(n_sub_n1.begin(),n_sub_n1.end(),J.begin(),J.end(),std::inserter(Jc, Jc.begin()));   
        // Lam2 = d2(J);
        for(int i = 0; i < J.size(); ++i){
            int temp = d2[J[i]];
            Lam2.push_back(temp);
        }
        // d3 = d2(Jc);
        vector<double>d3;
        for(int i = 0; i < Jc.size(); ++i){
            int temp = d2[Jc[i]];
            d3.push_back(temp);
        }
        // v3 = v2(Jc);
        vector<double>v3;
        for(int i = 0; i < Jc.size(); ++i){
            int temp = v2[Jc[i]];
            v3.push_back(temp);
        }
        
        int n2 = J.size();
        int n3 = n-n1-n2;
        // J = [J, Jc];
        for(int i = 0; i < Jc.size(); ++i){
            J.push_back(Jc[i]);
        }

        if(n3){
            rootfinder(d3,v3);
            /*
                [Lam3, tau, org, nflops1, percent] = rootfinder(d3, v3, N);
        
                [v3_hat, nflops1] = vhat(d3, Lam3, tau, org, v3, N);
            
                [s3, nflops1] = colnorms(d3, Lam3, tau, org, v3_hat, N);      
            */            
        }
    }
}