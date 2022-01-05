#include "superDC.h"
#include "BinTree.h"
#include "superdcmv_desc.h"
#include "superdcmv_node.h"
#include<string.h>
#include "eigenmatrix.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}
using namespace std;

void superdcmv_desc(EIG_MAT **Q,std::pair<int, int>*qSize, double **X,std::pair<int, int>xSize,BinTree *bt, int index, int ifTrans,std::pair<int,int>* indexRange,double N){
/*
%%% Input:
%%% Q: hss structured cauchylike eigenmatrix of descendants of i
%%% X: vectors to be multiplied
%%% ifTrans = 0: not transpose
%%% ifTranst = 1: transpose
%%% BinTree *bt: hss tree
%%% indexrange: indices range of each hss leaf block 


%%% Output
%%% ifTrans = 0: X = Q * X 
%%% ifTrans = 1: X = Q^T * X;
*/
    //smallest desendent of each node
    double * X_req = *X;
    std::vector<int> desc = bt->GetTreeDesc();  
    int k = desc[index+1];

    if(ifTrans == 0){
        //X = Q*X;
    }
    else{
        for(int j = k-1; j < index; ++j){
            
            int K = indexRange[j].second-indexRange[j].first+1; 
            int columns = xSize.second;
            double *tempX  = new double[K*columns];
            memset(tempX,0,sizeof(double)*K*columns);
            memcpy(tempX,X_req+indexRange[j].first,sizeof(double)*K*columns);
            std::pair<int,int>tempXSize = {K,columns};
            superdcmv_node(Q[j],qSize[j],&tempX,tempXSize,bt,j,1,N);
            memcpy(X_req+indexRange[j].first,tempX,sizeof(double)*K*columns);            
        }
    }
    *X = X_req;
    desc.clear();
    desc.shrink_to_fit();
}
