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
            if(A->bSizes[left-1].first <= A->bSizes[left-1].second){
                // D{c1} = D{c1} - norm(B{c1}) * (U{c1} * U{c1}');
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
                double B_c1_norm = norm_svd(tempB,A->bSizes[left-1]);
                
            }


        }

    }

}