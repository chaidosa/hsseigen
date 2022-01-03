#include "secular.h"
#include <assert.h>
#include <vector>
#include <set>
#include<algorithm>
#include<iostream>
#include <eigenmatrix.h>
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
    if(n1 < N){
       //conj(v2) ./ abs(v2); 
        
    }
    



}