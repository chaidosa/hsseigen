#include<bits/stdc++.h>
#include "../band2hss.h"
#include "../NPart.h"
#include "../BinTree.h"
using namespace std;

int main(){
    #ifndef INPUT_OUT
	    freopen("test_band_matrix.txt","r",stdin);	
    #endif
    long int n = 8192;
    int w = 9;
    long int m0 = 2048;
    double *A = new double[n*n];
    for(long int i = 0; i<n*n; i++){
        double temp;
        cin>>temp;
        A[i] = temp;
    }

    BinTree* bt=NULL;
	int *m=NULL;
	int mSize;
    int numNodes = 0;
	NPart(n, m0, &bt, &m, mSize, numNodes);

    B2HSS *res = new B2HSS();
    res = band2hss(&A, n*n, bt, m, mSize, w);

    cout<<"Generator creted successfully at:"<<res;
    double **temp;
    temp = res->D;
    delete[] temp;
    temp = res->U;
    delete[] temp;
    temp = res->R;
    delete[] temp;
    temp = res->B;
    delete[] temp;
    cout<<"Pritesh";
    delete[] A;
    return 0;
}