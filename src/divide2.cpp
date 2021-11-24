//divide2.cpp
//following code is stable version divide process in SuperDC.
#include<string.h>
#include<stack>
#include<vector>
#include <cmath>
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
    //smallest desendent of each node
    std::vector<int> desc = bt->GetTreeDesc();

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
        if(tch.size() == 0)
        {
            if(A->bSizes[left-1].first <= A->bSizes[left-1].second)
            {                
                double *tempD  = A->D[left-1];
                double *tempB  = A->B[left-1];
                double *tempU  = A->U[left-1];
                int size       = (A->uSizes[left-1].first)*(A->uSizes[left-1].second);
                double *tempC  = new double[size];

               //T = U{c1}*U{c1}'
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,A->uSizes[left-1].first,A->uSizes[left-1].second,A->uSizes[left-1].second,1,tempU,A->uSizes[left-1].second,tempU,A->uSizes[left-1].first,1,tempC,A->uSizes[left-1].second);

                double B_c1_norm = norm_svd(tempB,A->bSizes[left-1]);

                // D{c1} = D{c1} - norm(B{c1}) * (U{c1} * U{c1}');
                for(int row = 0 ; row < A->dSizes[left-1].first; row++)
                {
                    for(int col = 0; col < A->dSizes[left-1].second; col++)
                    {
                        A->D[left-1][col+row*(A->dSizes[left-1].second)] = A->D[left-1][col+row*(A->dSizes[left-1].second)] - B_c1_norm *(tempC[col+row*(A->dSizes[left-1].second)]);
                    }
                }
                delete [] tempC;
            }

            else
            {
                double *tempD  = A->D[left-1];
                double *tempB  = A->B[left-1];
                double *tempU  = A->U[left-1];
                double *tempt  = new double[A->uSizes[left-1].first * (A->bSizes[left-1].second)];
                double *temp   = new double[A->uSizes[left-1].first * (A->bSizes[left-1].second)];
                //T = U{c1} * B{c1};
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,A->uSizes[left-1].first,A->bSizes[left-1].second,A->uSizes[left-1].second,1,tempU,A->uSizes[left-1].second,tempB,A->bSizes[left-1].second,1,tempt,A->bSizes[left-1].second);
                //T*T'
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,A->uSizes[left-1].first,A->uSizes[left-1].first,A->uSizes[left-1].second,1,tempt,A->uSizes[left-1].first,tempt,A->uSizes[left-1].first,1,temp,A->uSizes[left-1].first);                
                double B_c1_norm = norm_svd(tempB,A->bSizes[left-1]);

                for(int row = 0 ; row < A->dSizes[left-1].first; row++)
                {
                    for(int col = 0; col < A->dSizes[left-1].second; col++)
                    {
                        A->D[left-1][col+row*(A->dSizes[left-1].second)] = A->D[left-1][col+row*(A->dSizes[left-1].second)] - temp[col+row*(A->dSizes[left-1].second)]/B_c1_norm;
                    }
                }

            }


        }
        else
        {
            if(A->bSizes[left-1].first <= A->bSizes[left-1].second)
            {
               std::vector<int> child = bt->GetChildren(left);           
               int eye_size           = A->rSizes[child[0]].second;
               double sqrt_B_c1_norm  = sqrt(norm_svd(A->B[left-1],A->bSizes[left-1]));
               double *Sp             = new double[eye_size*eye_size];
               memset(Sp,0,sizeof(*Sp)*eye_size*eye_size);
               for(int row_col = 0; row_col<eye_size; row_col++){
                   Sp[row_col+row_col*eye_size] = sqrt_B_c1_norm;
               } 

            }
            else
            {
                double *Sp             = new double[A->bSizes[left-1].first*A->bSizes[left-1].second];
                double sqrt_B_c1_norm  = sqrt(norm_svd(A->B[left-1],A->bSizes[left-1]));
                for(int row = 0; row < A->bSizes[left-1].first; row++)
                {
                    //memcpy(Sp+(row*A->bSizes[left-1].second), A->B[left-1]+row*A->bSizes[left-1].second, sizeof(double)*A->bSizes[left-1].second);
                    for(int col=0; col < A->bSizes[left-1].second; col++)
                    {
                        Sp[col+row*A->bSizes[left-1].second] = A->B[left-1][col+row*A->bSizes[left-1].second] / sqrt_B_c1_norm;  
                    }
                }
            }
            //Push(S, Sp), Push and Pop functionality yet to be implemented
            
            



        }

    }

}