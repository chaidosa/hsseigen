#include<cassert>
#include "bsxfun.h"
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


double* superdcmv(EIG_MAT **Qt, std::pair<int, int>*qSize, double *x, std::pair<int, int>* xSize, int * m, int mSize, int* I, BinTree* bt, int k, int ifTrans,double N){
	double* res=NULL;
    if(ifTrans == 1){
        assert(false);
    }
    else{
	    arrange_elements3(x, xSize->first, I, true); //x(I)=x if true. x=x(I) if false
            res = superdcmv_node(Qt[k],qSize[k],x,*xSize, bt,k,ifTrans,N);
    	    
	    int numNodes  = bt->GetNumNodes();
	    std::pair<int, int>* l = new std::pair<int,int>[numNodes];  
    	    for(int i = 0; i < numNodes; i++)
	    	l[i] = std::make_pair(0,0);    
    	
	    l[0] = {0,m[0] - 1};
    	    int it = 0;
    	    int lt = 0;

    	    for(int i = 0; i < numNodes; i++) {
		std::vector<int> ch = bt->GetChildren(i + 1);
		if(ch.size() == 0){
		    l[i] = {lt,lt + m[it] - 1};
		    lt   = l[i].second + 1;
		    it   = it + 1;
		}
		else{
		    l[i] = {l[ch[0] - 1].first,l[ch[1] - 1].second};
		}
	    }
    	    
	    res = superdcmv_desc(Qt,qSize,res,*xSize,bt,k,ifTrans,l,N);           
	    delete [] l;
    }
    return res;
}


