#include "superDC.h"
#include "BinTree.h"
#include "superdcmv_desc.h"
#include "superdcmv_node.h"
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
%%% ifTrans = 0: not transpose
%%% ifTranst = 1: transpose
%%% BinTree *bt: hss tree
%%% indexrange: indices range of each hss leaf block 


%%% Output
%%% ifTrans = 0: X = Q * X 
%%% ifTrans = 1: X = Q^T * X;
*/
    double alpha,beta;
    alpha = 1.0;
    beta  = 0.0;

    //smallest desendent of each node
    std::vector<int> desc = bt->GetTreeDesc();  
    int k = desc[index];

    if(ifTrans == 0){
        //X = Q*X;
    }
    else{
        for(int j = k-1; j <= index-1; ++j){
            double *req_X = X;
            int K = indexRange[j-1].second-indexRange[j-1].first+1; //Nikhil: Sure this is j-1? shouldn't this be j instead?
            int N = xSize.second;
            double *tempX  = new double[K*N];
            memset(tempX,0,sizeof(double)*K*N);
            //check
            double *tempX0;
            if(indexRange[j-1].first!=0) 
                tempX0=X+((indexRange[j-1].first+1)*N); //Nikhil: +1 is not necessary. You have stored X in row-major order and indexRange[x] gives you the pair [a,b] where a and b are the start and end indices of rows and columns of the square sub-matrix that node numbered j (on line 40) represents. So, row 0, e.g., starts at offset 0. Row 1 starts at offset N and so on..    
            else
                tempX0 = X;
	    //Nikhil: as you mentioned, this is not just cblas_dgemm for non-leaves. A series of updates to X is done (trying to find out the part in paper.) Normally, matrix vector product can be done with O(N^2) time. But structured matrix vector product is accelerated with the help of FMM modules here. superdcmv_node, superdcmv_cauchy, cauchylikematvec, and fmm1d_local_shift, fmm1d_local_shift_2 need to be implemented.
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,qSize[j-1].second,N,K,alpha,Q[j-1],qSize[j-1].second,tempX0,N,beta,tempX,N);
            memcpy(tempX0,tempX,sizeof(double)*qSize[j-1].second*N);
            delete [] tempX;
        }
    }
    desc.clear();
    desc.shrink_to_fit();
}
