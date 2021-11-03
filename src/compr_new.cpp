//Calculating QR using SVD
#include<string.h>
#include<algorithm>
#include<cmath>
#include<stdio.h>
#include<assert.h>
#include"compr.h"
#include"QR.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}

/*A is the input parameter initialized to T[i] for leaf i.
Q is the output parameter computing Ui, R is the output parameter computing Ri */
void compr_new(double* A, std::pair<int, int>aSize, double** Q, std::pair<int, int>& qSize, double** R, std::pair<int, int>& rSize, char* tol, double par)
{
    double *S    = new double[aSize.first*aSize.second];
    double *U    = new double[aSize.first*aSize.first];
    double *V    = new double[aSize.second*aSize.second];
    double *W    = new double[aSize.first*aSize.second];
    double *temp = new double[aSize.first*aSize.second];
    //calling lapack routine for SVD singular value decomposition
    LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'A','A',aSize.first,aSize.second,A,aSize.second,S,U,aSize.first,V,aSize.second,W);
    
    *Q = U;

    double *tempV = new double[aSize.first*aSize.second];

    //copying from V
    for(int i=0;i<aSize.first;i++){
        for(int j=0;j<aSize.second;j++){
            tempV[j+i*aSize.second]=V[j+i*aSize.second];
        }
    }

    //copying from S
    double *tempS = new double[aSize.first*aSize.first];

    for(int i=0;i<aSize.first;i++){
        for(int j=0;j<aSize.first;j++){
            tempS[j+i*aSize.first] = 0;
            if(i==j){
                tempS[j+i*aSize.first]=S[i];
            }
        }
    }

    double* r1 = S;
    double* r2 = V;
    delete [] r1;
    delete [] r2;
    S = tempS;
    V = tempV;

    //R = R*V'  Unlike Matlab, In this lapack routine V returned is already transposed so we don't need to transpose again.
    for(int i=0;i<aSize.first;i++)
    {
        for(int j=0;j<aSize.second;j++)
        {
            temp[j+i*aSize.second]=0;
            for(int k=0;k<aSize.first;k++)
            {
                temp[j+i*aSize.second]+=S[k+i*aSize.first]*(V[j+k*aSize.second]);
            }
        }
    }
    
    *R    = temp;
    qSize = {aSize.first,aSize.first};
    rSize = {aSize.first,aSize.second};
    delete [] S;
    delete [] W;
    delete [] V;

}