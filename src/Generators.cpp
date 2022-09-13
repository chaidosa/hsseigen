#include<string.h>
#include<math.h>
#include<algorithm>
#include<assert.h>
#include "BinTree.h"
#include "band2hss.h"
#include "Generators.h"
#include "iostream"
#include "test_mat2hsssym.h"

using namespace std;


GEN* HssGenerators(double *A, int aSize, BinTree* bt, int* m, int mSize, int w, int MorB){
    
    GEN *res;
    if(MorB == 2){
        cout << "Using Band2HSS\n";
        res = band2hss(A, aSize, bt, m, mSize, w);
    }
    else if(MorB == 1){
        cout << "Using Mat2Hssym\n";
        res = t_mat2hsssym(A, aSize, bt, m, mSize);      
    }
    else{
        cout<<"Something went wrong!!!, Try again\n";
        assert(false);
    }


    return res;
}