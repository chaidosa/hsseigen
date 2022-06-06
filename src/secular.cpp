#include "secular.h"
#include<bits/stdc++.h>
#include <assert.h>
#include <vector>
#include <set>
#include<algorithm>
#include<iostream>
#include "eigenmatrix.h"
#include "bsxfun.h"
#include "rootfinder.h"
#include "vhat.h"
#include  "colnorms.h"
using namespace std;

SECU* secular(double *d, int dSize, double *v, int vSize, double N){
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

    //Return values and global initializations
    nonleaf *res = new nonleaf();
    Root *rf_res=NULL; 
    SECU *ret = new SECU();
    //vector<double> v3_hat;
    double *v3_hat;
    int v3_hat_size;

    vector<double> tau;
    int *org;
    int org_size;
    double* s3;
    int s3_size=0;
    vector<double> d3;
    vector<double> v3;    
    vector<double> J,G;
    double *I;
    
    int n2;
    int n3;
    double percent;

    if(dSize != vSize){
        std::cout<<"diagonal and vector size should be same\n";
        assert(false);
    }
    //to store the results of secular 
    nonleaf *result = new nonleaf();

    double tol = 1.0e-6;
    int n      = vSize;
    //step 1: deflate small vi
    std::vector<double> Tempt;
    Tempt.reserve(vSize);
    std::vector<int> zero_2_n;
    zero_2_n.reserve(vSize);

    for(int i = 0; i < vSize; ++i)
    {

        if(std::abs(v[i]) < tol){
            Tempt.push_back(i);
        }
        //creating another array having indices 1:n
        zero_2_n.push_back(i);

    }

    std::vector<int>diff;
    diff.reserve(vSize);
    std::set_difference(zero_2_n.begin(),zero_2_n.end(),Tempt.begin(),Tempt.end(),std::inserter(diff, diff.begin()));

    std::vector<double> Lam1;
    Lam1.reserve(Tempt.size());

    for(int i = 0 ; i < Tempt.size() ; ++i){
        Lam1.push_back(d[(int)Tempt[i]]);
    }

    std::vector<double> Lam2;
    double* Lam3=NULL;

    std::vector<double>d2;
    d2.reserve(diff.size());
    std::vector<double>v2;
    v2.reserve(diff.size());
    int n1 = Tempt.size();

    for(int i = 0; i < diff.size(); ++i){
        d2.push_back(d[diff[i]]);
        v2.push_back(v[diff[i]]);    
        Tempt.push_back(diff[i]);
    }

// step 2: deflate  close eigenvalue
    double *v2c = new double[v2.size()];
    int v2c_size = v2.size();

    if(n1 < n){
       //conj(v2) ./ abs(v2);        
       for(int i = 0; i < v2.size(); ++i){
           v2c[i] = v2[i] / abs(v2[i]);
           v2[i] = abs(v2[i]);
       }
       
       //sort array by keeping track of previous index
       
       //vector to store element       
        vector<pair<double,int>> vP;
        vP.reserve(d2.size());
        for(int i = 0; i < d2.size(); ++i){
            vP.push_back(make_pair(d2[i],i));
        }

        //sorting pair vector
        sort(vP.begin(),vP.end());

        //to store the new sorted indices
        I = new double[d2.size()];

        //copying sorted values back to d2;
        for(int i = 0; i < d2.size(); i++){
            d2[i] = vP[i].first;
            I[i] = vP[i].second;
        }

        vector<double> tempV2(v2);
        for(int i = 0 ;i < v2.size(); ++i)
            v2[i] = tempV2[vP[i].second];        

        int p=0;        
        for(int j = 1; j < d2.size(); j++){
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
   
        
        vector<int> n_sub_n1;
        n_sub_n1.reserve(n-n1);
        for(int i = 0; i < (n-n1); ++i){
            n_sub_n1.push_back(i);
        }
        
        std::vector<int>Jc;
        Jc.reserve(n-n1);
        std::set_difference(n_sub_n1.begin(),n_sub_n1.end(),J.begin(),J.end(),std::inserter(Jc, Jc.begin())); 

        // Lam2 = d2(J);
        for(int i = 0; i < J.size(); ++i){
            double temp = d2[J[i]];
            Lam2.push_back(temp);
        }

        // d3 = d2(Jc); 
        
        d3.reserve(Jc.size());
        v3.reserve(Jc.size());   
        
        n2 = J.size();
        n3 = n-n1-n2;
        
        for(int i = 0; i < Jc.size(); ++i){
            double temp = d2[Jc[i]];
            d3.push_back(temp);
            
            temp = v2[Jc[i]];
            v3.push_back(temp);
            J.push_back(Jc[i]);
        }
        // J = [J, Jc];

        if(n3){
            //Rootfinder call            
            rf_res = rootfinder(d3, v3, N);
            Lam3 = rf_res->x;
            tau = rf_res->tau;
            org = rf_res->org;
            org_size = rf_res->org_size;
            percent = rf_res->percent;
            //V hat call
            v3_hat = vhat(d3, Lam3, org, org_size, tau, v3);

            //Colnorms call
            s3 = colnorms(d3, Lam3, tau, org, org_size, v3_hat); 
	        s3_size=org_size;
        }
        
    }

    int lamSize = Lam1.size() + Lam2.size() + org_size;
    double *Lam = new double[lamSize];   
   

    double *pointLam = Lam;
    memcpy(pointLam, Lam1.data(), sizeof(double)*Lam1.size());
    pointLam = Lam + Lam1.size();
    memcpy(pointLam, Lam2.data(), sizeof(double)*Lam2.size());
    pointLam = Lam + Lam1.size() + Lam2.size();
    if(Lam3)
    memcpy(pointLam, Lam3, sizeof(double)*org_size); 

    if(n1 < n){    
    
        
        res->QC[0] = v3_hat;
        res->qcSizes[0] = make_pair(org_size, 1);

	    res->QC[1] = s3;
        res->qcSizes[1] = make_pair(s3_size, 1);


        double *d33 = new double[d3.size()];
        memcpy(d33, &d3[0], sizeof(double)*d3.size());
        res->QC[2] = d33;
        res->qcSizes[2] = {d3.size(), 1};


        /*double *Lam33 = new double[Lam3.size()];
        memcpy(Lam33, &Lam3[0], sizeof(double)*Lam3.size());
        res->QC[3] = Lam33;
        res->qcSizes[3] = {Lam3.size(), 1};*/
	    res->QC[3]=Lam3;
	    if(Lam3)
        	res->qcSizes[3] = std::make_pair(org_size,1);
	    else
        	res->qcSizes[3] = std::make_pair(0,1);



        double *tau1 = new double[tau.size()];
        memcpy(tau1, &tau[0], sizeof(double)*tau.size());
        res->QC[4] = tau1;
        res->qcSizes[4] = {tau.size(), 1};


        //double *org1 = new double[org.size()];
        //memcpy(org1, org.data(), sizeof(double)*org.size());
        //res->QC[5] = org1;
        res->Org = org;
        res->qcSizes[5] = {org_size, 1};

        res->I = I;
        res->ISize = {d2.size(),1};
    
        double *G1 = new double[G.size()];
        memcpy(G1, &G[0], sizeof(double)*G.size());
        res->G = G1;    
        res->GSize = {G.size(), 1};

        double *J1 = new double[J.size()];
        memcpy(J1, &J[0], sizeof(double)*J.size());
        res->J = J1;
        res->JSize = {J.size(), 1};

        double *T = new double[Tempt.size()];
        memcpy(T, Tempt.data(), sizeof(double)*Tempt.size());
        res->T = T;
        res->TSize = {Tempt.size(), 1};

        res->v2c = v2c;
        res->v2cSize = {v2c_size,1};   


        res->n = n;
        res->n1 = n1;
        res->n2 = n2;
        res->n3 = n3;
        ret->Q = res;
        ret->Lam = Lam;
    
        return ret;
    }

    else
    {        
        if(Lam3)
            delete [] Lam3;
        res->JSize = {0,0};
        res->GSize = {0,0};
        res->ISize = {0,0};
        
        for(int i = 0; i < 6; i++)
            res->qcSizes[i] = make_pair(0,0);
        
        res->v2cSize = make_pair(0,0);
        double *T = new double[Tempt.size()];
        memcpy(T, &Tempt[0], sizeof(double)*Tempt.size());
        
        res->T = T;
        res->TSize = {Tempt.size(), 1};
        res->n = n;
        res->n1 = n1;
        res->n2 = 0;
        res->n3 = 0;
        ret->Q = res;
        ret->Lam = Lam;
        ret->percent = percent;
    }

    return ret;
}
