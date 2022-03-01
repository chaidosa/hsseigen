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
    std::pair<int, int>rg[bt->GetNumNodes()];

    for(int itr = 0; itr < bt->GetNumNodes(); itr++){
        rg[itr]  = {indexRange[itr].first-indexRange[index].first, indexRange[itr].second-indexRange[index].first};
    }

    if(ifTrans == 0){
        //X = Q*X;
    }
    
    else{
        for(int j = k-1; j < index; ++j){
            
            int K = rg[j].second-rg[j].first+1; 
            int columns = xSize.second;
            double *tempX  = new double[K*columns];
            //memset(tempX,0,sizeof(double)*K*columns);
            memcpy(tempX,X_req+(rg[j].first*columns),sizeof(double)*K*columns);
            
            std::pair<int,int>tempXSize = {K,columns};

            superdcmv_node(&Q[j],qSize[j],&tempX,tempXSize,bt,j,1,N);

            memcpy(X_req+(rg[j].first * columns),tempX,sizeof(double)*K*columns);     
            delete [] tempX;       
        }
    }
    *X = X_req;
    desc.clear();
    desc.shrink_to_fit();
}
