#include "superDC.h"
#include "BinTree.h"
#include "superdcmv_desc.h"
#include<string.h>
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}
using namespace std;

void superdcmv_desc(double **Q,std::pair<int, int>*qSize, double *X,std::pair<int, int>xSize,BinTree *bt, int index, int ifTrans,std::pair<int,int>* indexRange){
/*
%%% Input:
%%% Q: hss structured cauchylike eigenmatrix of descendants of i
%%% X: vectors to be multiplied
%%% t = 0: not transpose
%%% t = 1: transpose
%%% BinTree *bt: hss tree
%%% rg: indices range of each hss leaf block 
%%% desc: smallest descendant of node i in the hss tree tr

%%% Output
%%% t = 0: X = Q * X 
%%% t = 1: X = Q^T * X;
*/
    double alpha,beta;
    alpha = 1.0;
    beta  = 0.0;
    std::vector<int> desc = bt->GetTreeDesc();  
    int k = desc[index];

    if(ifTrans == 0){
        //X = Q*X;
    }
    else{
        for(int j = k; k <= index-1; ++k){
            double *req_X = X;
            int K = indexRange[j-1].second-indexRange[j-1].first+1;
            int N = xSize.second;
            double *tempX  = new double[K*N];
            memset(tempX,0,sizeof(double)*K*N);

            double *tempX0;
            if(indexRange[j-1].first!=0)
                tempX0=X+((indexRange[j-1].first+1)*N);
            else
                tempX0 = X;

            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,qSize[j-1].second,N,K,alpha,Q[j-1],qSize[j-1].second,tempX0,N,beta,tempX,N);
            memcpy(tempX0,tempX,sizeof(double)*qSize[j-1].second*N);
            delete [] tempX;
        }
    }
    desc.clear();
    desc.shrink_to_fit();
}
