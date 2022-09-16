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


double* superdcmv(EIG_MAT **Qt, std::pair<int, int>*qSize, double *x, std::pair<int, int> xSize, BinTree* bt, int* m, int ifTrans,double N){
	double* res=NULL;
    if(ifTrans == 0){
        assert(false);
    }
    else{
    	    int K = bt->GetNumNodes();
	    std::pair<int, int>* l = new std::pair<int,int>[K];  

	    for(int k = 0; k < K; k++)
		l[k] = std::make_pair(0,0);    

	    l[0]                  = {0,m[0] - 1};
	    int it                = 0;
	    int lt                = 0;

	    for(int k = 0; k < K; k++)
	    {
		std::vector<int> ch = bt->GetChildren(k + 1);
		if(ch.size() == 0)
		{
		    l[k] = {lt,lt + m[it] - 1};
		    lt   = l[k].second + 1;
		    it   = it + 1;
		}
		else
		    l[k] = {l[ch[0] - 1].first,l[ch[1] - 1].second};
	    }


	    std::pair<int, int> qSizek=qSize[K-1];
            res = superdcmv_node(Qt[K-1],xSize,x,xSize, bt,K-1,1,N);
    	    res = superdcmv_desc(Qt,qSize,res,xSize,bt,K-1,1,l,N);           
        }
    return res;
}


