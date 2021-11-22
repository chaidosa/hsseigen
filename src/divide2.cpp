//divide2.cpp :: Pritesh Verma
//following code is stable version divide process in SuperDC.

#include "test_mat2hsssym.h"
#include "BinTree.h"
#include "divide2.h"
#include<stdio.h>
#include "QR.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}
double norm_svd(double * A, std::pair<int, int>aSize){

    double *S = new double[aSize.first*aSize.first];
    double *U,*V,*su;
    LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'N','N',aSize.first,aSize.second,A,aSize.second,S,U,aSize.first,V,aSize.second,su);
    return S[0];
}


DVD* divide2(tHSSMat *A, BinTree *bt,int* m, int mSize){

    int n                 = bt->GetNumNodes();

    std::vector<int> desc = bt->GetDescendents(n);

    //Dividing process
    for(int i = n;i >= 1; i--){
    
        std::vector<int> ch = bt->GetChildren(i);

        if(ch.size() == 0)
            continue;
        
        //divide at the non-leaf nodes
        int left  = ch[0];
        int right = ch[1];

        //updating B generators of left subtree
        std::vector<int> tch = bt->GetChildren(left);
        if(tch.size() == 0){
            if(A->bSizes[left-1].first <= A->bSizes[left-1].second)
            {                
                double *tempD  = A->D[left-1];
                double *tempB  = A->B[left-1];
                double *tempU  = A->U[left-1];
                int size       = (A->uSizes[left-1].first)*(A->uSizes[left-1].second);
                double *tempUt = new double[size];
                double *tempC  = new double[size];
                //use memcopy here
                for(int i = 0; i < size; i++)
                    tempUt[i] = tempU[i];
                
                GetTransposeInPlace(tempUt, A->uSizes[left-1].first, A->uSizes[left-1].second);
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,A->uSizes[left-1].first,A->uSizes[left-1].second,A->uSizes[left-1].second,1,tempU,A->uSizes[left-1].second,tempUt,A->uSizes[left-1].first,1,tempC,A->uSizes[left-1].second);
                //norm(B{c1})
                double B_c1_norm = norm_svd(tempB,A->bSizes[left-1]);

                // D{c1} = D{c1} - norm(B{c1}) * (U{c1} * U{c1}');
                for(int row = 0 ; row < A->dSizes[left-1].first; row++){
                    for(int col = 0; col < A->dSizes[left-1].second; col++){
                        A->D[left-1][col+row*(A->dSizes[left-1].second)] = A->D[left-1][col+row*(A->dSizes[left-1].second)] - B_c1_norm *(tempC[col+row*(A->dSizes[left-1].second)]);
                    }
                }
                delete [] tempC;
                delete [] tempUt;
            }

            else
            {
                double *tempD  = A->D[left-1];
                double *tempB  = A->B[left-1];
                double *tempU  = A->U[left-1];
                double *tempt  = new double[A->uSizes[left-1].first * (A->bSizes[left-1].second)];

                //T = U{c1} * B{c1};
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,A->uSizes[left-1].first,A->bSizes[left-1].second,A->uSizes[left-1].second,1,tempU,A->uSizes[left-1].second,tempB,A->bSizes[left-1].second,1,tempt,A->bSizes[left-1].second);
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,A->uSizes[left-1].first,A->uSizes[left-1].first,A->uSizes[left-1].second,1,tempt,A->uSizes[left-1].first,tempt,A->uSizes[left-1].first,1,tempt,A->uSizes[left-1].first);
                
            }


        }

    }

}