#include<math.h>
#include<assert.h>
#include"BinTree.h"


//This function creates a perfect binary tree with n nodes.
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



//This function creates a full binary tree with n nodes.
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



//This function creates full binary tree with nodes n.
void BinTree::Create(int n){
    
    if(n%2 == 0){
        std::cout<<"Error!! cannot create the tree input must be odd";
        assert(false);
    }

    numNodes = n;
    std::cout<<"Creating tree with "<<n<<" nodes";

    tr = ntree(n);
    for(int i = 0; i < numNodes; i++){
        if(tr[i] == 0)
            continue;
        ch[tr[i]].push_back(i+1);    
    }

    nodeAtLvl = hsslevel();
    //ntree(int n);
    
}



//This function returns the childrens of a node. 
std::vector<int> BinTree::GetChildren(int ID)
{
	std::vector<int> ret;
    std::map<int , std::vector<int>>::iterator it;
    it = ch.find(ID);
    if(it != ch.end())
        ret.insert(ret.begin(), it->second.begin(), it->second.end());

	return ret;
}



//This function returns the smallest descendents of all the nodes.
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



// This function returns number of nodes in the tree.
int BinTree::GetNumNodes(){
    return numNodes;
}



//This function returns all the nodes stored in level order.
std::vector<vector<int>> BinTree::hsslevel(){
    
    std::vector<vector<int>>nodes_lvl;
    
    std::vector<int> curr;
    curr.push_back(numNodes);
    while(!curr.empty())
    {        
        nodes_lvl.push_back(curr);

        vector<int>nextlvl;

        for(int j = 0; j < curr.size(); j++){
            vector<int>ch = GetChildren(curr[j]);
            if(!ch.empty())
                nextlvl.insert(nextlvl.end(), ch.begin(), ch.end());
        }

        curr.clear();
        curr.insert(curr.begin(), nextlvl.begin(), nextlvl.end());
        nextlvl.clear();
    }  
    numLevels = nodes_lvl.size();
    return nodes_lvl;
}
