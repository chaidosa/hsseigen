#include "BinTree.h"
#include "test_mat2hsssym.h"
#include "superDC.h"
#include "divide2.h"
SDC* superDC(tHSSMat* A,  BinTree* bt, int* m, int mSize){

//Dividing Stage
DVD *resDvd = divide2(A,bt,m,mSize);


int N  = bt->GetNumNodes();
//Index range of each node
std::pair<int,int>* l = new std::pair<int,int>[N];

    for(int i = 0; i < N; i++)
		l[i] = std::make_pair(0,0);    

    l[0]                  = {0,m[0]-1};
    int it                = 0;
    int lt                = 0;

    for(int i = 0; i < N; i++)
    {
        std::vector<int> ch = bt->GetChildren(i+1);
        if(ch.size() == 0)
        {
            l[i] = {lt,lt+m[it]-1};
            lt   = l[i].second +1;
            it   = it+1;
        }
        else
        {
            l[i] = {l[ch[0]-1].first,l[ch[1]-1].second};

        }
    }    

    


} 
