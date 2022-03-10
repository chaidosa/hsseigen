#include<math.h>
#include<assert.h>
#include"BinTree.h"

std::vector<int> btree(int n){
    
    if(n%2 == 0){
        std::cout<<"Error!! cannot create the tree input must be odd";
        assert(false);
    }

    if(n > 3){
        int m = (n-1)/2;
        std::vector<int> t11 = btree(m);
        std::vector<int> t22(t11.begin(), t11.end());
        
        for(int i = 0; i < t22.size(); i++)
            t22[i] = t22[i] + m;

        t11[t11.size() - 1] = n;
        t22[t22.size() - 1] = n;
        
        std::vector<int> result(t11);
        result.insert(result.end(), t22.begin(), t22.end());      

        result.push_back(0);

        return result;
    }

    else if (n == 3){
        std::vector<int>result;
        result.push_back(3);
        result.push_back(3);
        result.push_back(0);
        return result;
    }
    
    else{
        std::vector<int> result;
        result.push_back(0);
        return result;
    }

}


std::vector<int> ntree(int n){
    
    if(n%2 == 0){
        std::cout<<"Error!! cannot create the tree input must be odd";
        assert(false);
    }

    if(n == 1){
        std::vector<int> result;
        result.push_back(0);
        return result;
    }    

    else{
        
        int n1 = std::floor(std::log2(n));
            n1 = std::pow(2, n1) - 1;

        std::vector<int>tr1 = btree(n1); 

        std::vector<int>tr2 = ntree(n-n1-1);
        for(int i = 0; i < tr2.size(); i++)
            tr2[i] = tr2[i] + n1;
        
        tr1[tr1.size() - 1] = n;
        tr2[tr2.size() - 1] = n;      


        std::vector<int>result(tr1);
        result.insert(result.end(), tr2.begin(), tr2.end());
        result.push_back(0);
        return result;
    }



}
void BinTree::Create(int n){
    
    if(n%2 == 0){
        std::cout<<"Error!! cannot create the tree input must be odd";
        assert(false);
    }

    numNodes = n;
    std::cout<<"Creating tree with "<<n<<" nodes";

    tr = ntree(n);   
    //ntree(int n);
    
}

std::vector<int> BinTree::GetChildren(int ID)
{
	std::vector<int> ret;
	for(int i=0;i<numNodes;i++)
	{
		if(tr[i] == ID)
		{
			ret.push_back(i+1);
			if(ret.size() == 2)
				break;
		}
	}

	return ret;
}

std::vector<int> BinTree::GetTreeDesc(){
    std::vector<int> result(numNodes+1, 0);

    for(int i = 1; i < numNodes+1; i++){
        if(GetChildren(i).size() == 0){
            result[i] = i;
        }
        else{
            result[i] = result[GetChildren(i)[0]];
        }
    }

    return result;
}


int BinTree::GetNumNodes(){
    return numNodes;
}
