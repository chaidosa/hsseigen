#include<cassert>
#include "superdcmv.h"
#include "superdcmv_node.h"
#include "superdcmv_desc.h"

/*[Q0, I, tr, m] = Q{:};
k = length(tr);
rg = indrange(tr, m);
desc = treedesc(tr);
nflops = 0;

if t == 0
    X(I, :) = X;
    [X, nflops1] = superdcmv_node(Q0{k}, X, tr, k, t, N);
    [X, nflops2] = superdcmv_desc(Q0, X, tr, k, t, rg, desc, N);
    nflops = nflops + nflops1 + nflops2;
else
    [X, nflops1] = superdcmv_desc(Q0, X, tr, k, t, rg, desc, N);
    [X, nflops2] = superdcmv_node(Q0{k}, X, tr, k, t, N);
    X = X(I, :);
    nflops = nflops + nflops1 + nflops2;
end*/


double* superdcmv(EIG_MAT **Qt, std::pair<int, int>*qSize, double *x, BinTree* bt, int ifTrans,double N){
	double* res=NULL;
    if(ifTrans == 0){
        assert(false);
    }
    else{
            int columns = qSize->second;
            std::pair<int,int>tempXSize = {columns,1};
    	    int k = bt->GetNumNodes();
            res = superdcmv_node(Qt[k],qSize[k],x,tempXSize, bt,k,1,N);
	    std::pair<int, int> l={0, columns-1};
    	    res = superdcmv_desc(Qt,qSize,res,tempXSize,bt,k,1,&l,N);           
        }
    return res;
}


