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
void compr_new(double* A, std::pair<int, int>aSize, double** Q, std::pair<int, int>& qSize, double** R, std::pair<int, int>& rSize, char const* tol, double par)
{
    int rk = std::min(aSize.first,aSize.second);
    double *S    = new double[rk];
    double *U    = new double[aSize.first*(aSize.first)];
    double *V    = new double[aSize.second*(aSize.second)];
    double *W    = new double[aSize.first*(aSize.second)];
   
    //calling lapack routine for SVD singular value decomposition
    LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'A','A',aSize.first,aSize.second,A,aSize.second,S,U,aSize.first,V,aSize.second,W);    
   
    double *tempU = new double[aSize.first*rk];
    //for(int row = 0; row<aSize.first; row++)
    //    memcpy(tempU+row*rk,U+row*aSize.first,sizeof(double)*rk);
    for(int row =0 ;row<aSize.first; row++){
        for(int col=0; col<rk; col++){
            tempU[col+(row*rk)] = U[col+(row*(aSize.first))];
        }
    }

    *Q = tempU;
    tempU = NULL;
    qSize = {aSize.first,rk};

    double *tempV = new double[rk*(aSize.second)];
    double *temp = new double[rk*(aSize.second)];
    //copying from V
    for(int i=0;i<rk;i++){
        for(int j=0;j<aSize.second;j++){
            tempV[j+(i*(aSize.second))]=V[j+i*(aSize.second)];
        }
    }

    //copying from S
    double *tempS = new double[rk*rk];

    for(int i=0;i<rk;i++){
        for(int j=0;j<rk;j++){
            tempS[j+i*rk] = 0;
            if(i==j){
                tempS[j+i*rk]=S[i];
            }
        }
    }

 //  double* r1 = S;
 //   double* r2 = V;
    
 //  S = tempS;
  //  V = tempV;

    //R = R*V'  Unlike Matlab, In this lapack routine V returned is already transposed so we don't need to transpose again.
    for(int i=0;i<rk;i++)
    {
        for(int j=0;j<aSize.second;j++)
        {
            temp[j+i*aSize.second]=0;
            for(int k=0;k<rk;k++)
            {
                temp[j+i*(aSize.second)]+=tempS[k+i*rk]*(tempV[j+k*(aSize.second)]);
            }
        }
    }
    
    *R    = temp;
    temp  = NULL;   
    rSize = {rk,aSize.second};
    delete [] S;
    delete [] W;
    delete [] V;
    delete [] U;
    delete [] tempS;
    delete [] tempV;
  //  delete [] r1;
  //  delete [] r2;
    
  

}
