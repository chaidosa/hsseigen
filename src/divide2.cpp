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
#include "band2hss.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}
using namespace std;
//norm using svd
void norm_svd(double * A, std::pair<int, int>aSize, double *norm){

    double *S  = new double[aSize.first*aSize.second];
    double *su = new double[aSize.first*aSize.second];
    double *U  = new double[aSize.first*aSize.first];
    double *V  = new double[aSize.second*aSize.second];
    int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'N','N',aSize.first,aSize.second,A,aSize.second,S,U,aSize.first,V,aSize.second,su);
    if( info > 0 ) {
                printf("ERROR: LAPACKE_dgesvd failed to converge.\n" );
                exit( 1 );
    }
    delete [] U;
    delete [] V;
    delete [] su;
    *norm = S[0];
    delete [] S;
}


DVD* divide2(tHSSMat *A, BinTree *bt,int* m, int mSize)
{
    DVD* ret = new DVD();
    double alpha,beta;
    alpha = 1.0;
    beta  = 0.0;
    int n                 = bt->GetNumNodes();
    //smallest desendent of each node
    std::vector<int> desc = bt->GetTreeDesc();   

    //Dividing process
    for(int i = n;i >= 1; i--)
    {
    
        std::vector<int> ch = bt->GetChildren(i);

        if(ch.empty())
            continue;
        
        //divide at the non-leaf nodes
        int left  = ch[0];
        int right = ch[1];

        double sqrt_B_c1_norm=0;
        double B_c1_norm = 0;
        double *Temp_BC1 = new double[(A->bSizes[left-1].first)*(A->bSizes[left-1].second)];
        memcpy(Temp_BC1,A->B[left-1],sizeof(double)*(A->bSizes[left-1].first)*(A->bSizes[left-1].second));
        norm_svd(Temp_BC1,A->bSizes[left-1],&B_c1_norm);
        sqrt_B_c1_norm=std::sqrt(B_c1_norm);
        delete [] Temp_BC1;
        //updating B generators of left subtree
        std::vector<int> tch = bt->GetChildren(left);
        if(tch.empty())
        {
            if(A->bSizes[left-1].first <= A->bSizes[left-1].second)
            {                
                double *tempU  = A->U[left-1];                
                double *tempUUt  = new double[A->uSizes[left-1].first*(A->uSizes[left-1].first)];
                memset(tempUUt,0,sizeof(double)*(A->uSizes[left-1].first*(A->uSizes[left-1].first)));
               //T = U{c1}*U{c1}'
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,A->uSizes[left-1].first,A->uSizes[left-1].first,A->uSizes[left-1].second,alpha,tempU,A->uSizes[left-1].second,tempU,A->uSizes[left-1].second,beta,tempUUt,A->uSizes[left-1].first);
                               
                // D{c1} = D{c1} - norm(B{c1}) * (U{c1} * U{c1}');
                for(int row = 0 ; row < A->dSizes[left-1].first; row++)
                {
                    for(int col = 0; col < A->dSizes[left-1].second; col++)
                    {
                        A->D[left-1][col+row*(A->dSizes[left-1].second)] = A->D[left-1][col+row*(A->dSizes[left-1].second)] - (B_c1_norm *(tempUUt[col+row*(A->uSizes[left-1].first)]));
                    }
                }
                delete [] tempUUt;
            }

            else
            {
                double *tempB  = A->B[(left-1)];
                double *tempU  = A->U[(left-1)];
                double *tempT  = new double[(A->uSizes[left-1].first)*(A->bSizes[left-1].second)];
                double *tempTTt   = new double[(A->uSizes[left-1].first)*(A->uSizes[left-1].first)];
                memset(tempT,0,sizeof(double)*(A->uSizes[left-1].first)*(A->bSizes[left-1].second));
                memset(tempTTt,0,sizeof(double)*(A->uSizes[left-1].first)*(A->uSizes[left-1].first));
                //T = U{c1} * B{c1};
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,A->uSizes[left-1].first,A->bSizes[left-1].second,A->uSizes[left-1].second,alpha,tempU,A->uSizes[left-1].second,tempB,A->bSizes[left-1].second,beta,tempT,A->bSizes[left-1].second);
                //T*T'
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,A->uSizes[left-1].first,A->uSizes[left-1].first,A->bSizes[left-1].second,alpha,tempT,A->bSizes[left-1].second,tempT,A->bSizes[left-1].second,beta,tempTTt,A->uSizes[left-1].first);                
                
                for(int row = 0 ; row < A->dSizes[left-1].first; row++)
                {
                    for(int col = 0; col < A->dSizes[left-1].second; col++)
                    {
                        A->D[left-1][col+row*(A->dSizes[left-1].second)] = A->D[left-1][col+row*(A->dSizes[left-1].second)] - (tempTTt[col+row*(A->uSizes[left-1].first)]/B_c1_norm);
                    }
                }
                delete [] tempTTt;
                delete [] tempT;
            }
        } 

        else
        {
            //temporary vector to push into the stack            
            double *Sp;
            std::pair<int,int>tempP;
            if(A->bSizes[left-1].first <= A->bSizes[left-1].second)
            {
                std::vector<int> child = bt->GetChildren(left);           
                int eye_size           = A->rSizes[(child[0]-1)].second;               
                Sp                    = new double[eye_size*eye_size];
                memset(Sp,0,sizeof(*Sp)*eye_size*eye_size);

                for(int row_col = 0; row_col<eye_size; row_col++)
                {
                   Sp[row_col+row_col*eye_size] = sqrt_B_c1_norm;
                }             
                
                tempP = {eye_size,eye_size};

            }

            else
            {
                Sp = new double[(A->bSizes[left-1].first)*(A->bSizes[left-1].second)];
                memset(Sp,0,sizeof(double)*(A->bSizes[left-1].first)*(A->bSizes[left-1].second));
               
                for(int row = 0; row < A->bSizes[left-1].first; row++)
                {
                    for(int col=0; col < A->bSizes[left-1].second; col++)
                    {
                        Sp[col+row*(A->bSizes[left-1].second)] = A->B[left-1][col+row*(A->bSizes[left-1].second)] / sqrt_B_c1_norm;  
                    }
                }

                tempP = {(A->bSizes[left-1].first),(A->bSizes[left-1].second)};
            }
           
            std:: stack<pair<double*,pair<int,int>>> S;
            S.push({Sp,tempP});

            for(int j=left-1;j >= desc[left];j--)
            {
                Sp         = S.top().first;               
                int Sp_row = (S.top().second).first;
                int Sp_col = (S.top().second).second;
                S.pop();

                std::vector<int> temp_ch = bt->GetChildren(bt->tr[j-1]);
                if(j== temp_ch[1])
                {
                    int sib = temp_ch[0];                    
                    S.push({Sp,{Sp_row,Sp_col}});

                    //Sp*Sp'
                    double *tempSST = new double[Sp_row*Sp_row];
                    memset(tempSST,0,sizeof(double)*Sp_row*Sp_row);                    
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,Sp_row,Sp_row,Sp_col,alpha,Sp,Sp_col,Sp,Sp_col,beta,tempSST,Sp_row);

                    //R{sib}*(Sp*Sp')
                    double *R_sib = A->R[(sib-1)]; 
                    double *tempRSST = new double[(A->rSizes[(sib-1)].first)*Sp_row];
                    memset(tempRSST,0,sizeof(double)*(A->rSizes[(sib-1)].first)*Sp_row);
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(A->rSizes[(sib-1)].first),Sp_row,Sp_row,alpha,R_sib,(A->rSizes[(sib-1)].second),tempSST,Sp_row,beta,tempRSST,Sp_row);

                    double *tempRSST_Rjt = new double[(A->rSizes[(sib-1)].first)*(A->rSizes[(j-1)].first)];
                    memset(tempRSST_Rjt,0,sizeof(double)*(A->rSizes[(sib-1)].first)*(A->rSizes[(j-1)].first));
                    double *Rj = A->R[j-1];
                    
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(A->rSizes[(sib-1)].first),(A->rSizes[(j-1)].first),(A->rSizes[(j-1)].second),alpha,tempRSST,Sp_row,Rj,(A->rSizes[(j-1)].second),beta,tempRSST_Rjt,(A->rSizes[(j-1)].first));

                    for(int row = 0; row<A->bSizes[sib-1].first;row++)
                    {
                        for(int col=0; col<A->bSizes[sib-1].second;col++)
                        {
                            A->B[sib-1][col+(row*(A->bSizes[sib-1].second))] = A->B[sib-1][col+(row*(A->bSizes[sib-1].second))] - tempRSST_Rjt[col+row*A->rSizes[j-1].first];
                        }
                    }
                    delete [] tempSST;
                    delete [] tempRSST;
                    delete [] tempRSST_Rjt;
                }
                //Sj = R{j} * Sp;
                double *Sj = new double[A->rSizes[j-1].first*Sp_col];
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,A->rSizes[j-1].first,Sp_col,Sp_row,alpha,A->R[j-1],A->rSizes[j-1].second,Sp,Sp_col,beta,Sj,Sp_col); 

                std::vector<int> temp_ch2 = bt->GetChildren(j);
                if(temp_ch2.size()==0)
                {
                    //T = U{j}*Sj;
                    double *tempT = new double[A->uSizes[j-1].first*Sp_col];
                    memset(tempT,0,sizeof(double)*(A->uSizes[j-1].first)*Sp_col);                    
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,A->uSizes[j-1].first,Sp_col,A->uSizes[j-1].second,alpha,A->U[j-1],(A->uSizes[j-1].second),Sj,Sp_col,beta,tempT,Sp_col);
                    //T*T'
                    double *tempTTt = new double[(A->uSizes[j-1].first)*(A->uSizes[j-1].first)]; 
                    memset(tempTTt,0,sizeof(double)*(A->uSizes[j-1].first)*(A->uSizes[j-1].first));  
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,A->uSizes[j-1].first,A->uSizes[j-1].first,Sp_col,alpha,tempT,Sp_col,tempT,Sp_col,beta,tempTTt,A->uSizes[j-1].first);
                    // D{j}  = D{j} - T * T';
                    for(int row =0;row < (A->dSizes[j-1].first);row++)
                    {
                        for(int col=0;col < (A->dSizes[j-1].second);col++)
                        {
                            A->D[j-1][col+(row*(A->dSizes[j-1].second))] -= tempTTt[col+(row*(A->uSizes[j-1].first))];
                        }
                    }
                    delete [] tempT;
                    delete [] tempTTt;
                    delete [] Sj;
                }
                else
                {
                    S.push({Sj,{A->rSizes[j-1].first,Sp_col}});
                }
               // delete[] Sj;
               // delete[] Sp;
            }
            //If stack is not empty
            while (!S.empty())
            {
                double *temp = S.top().first;
                S.pop();
                delete [] temp;
                temp=NULL;
            }          

        } //end of update of left subtree

        //Update B generators of right subtree
        tch = bt->GetChildren(right);
        if(tch.size()==0){
            if(A->bSizes[left-1].first <= A->bSizes[left-1].second)
            {
                //U{c2}*B{c1}'
                double *tempT = new double[(A->uSizes[right-1].first)*(A->bSizes[left-1].first)];
                memset(tempT,0,sizeof(double)*(A->uSizes[right-1].first)*(A->bSizes[left-1].first));
                double *tempU = A->U[right-1];
                double *tempB = A->B[left-1];
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(A->uSizes[right-1].first),(A->bSizes[left-1].first),(A->bSizes[left-1].second),alpha,tempU,(A->uSizes[right-1].second),tempB,(A->bSizes[left-1].second),beta,tempT,(A->bSizes[left-1].first));
                //T*T'
                double *tempTTt = new double[(A->uSizes[right-1].first)*(A->uSizes[right-1].first)];
                memset(tempTTt,0,sizeof(double)*(A->uSizes[right-1].first)*(A->uSizes[right-1].first));
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(A->uSizes[right-1].first),(A->uSizes[right-1].first),(A->bSizes[left-1].first),alpha,tempT,(A->bSizes[left-1].first),tempT,(A->bSizes[left-1].first),beta,tempTTt,(A->uSizes[right-1].first));
                              
                // D{c2} = D{c2} - (T * T') / norm(B{c1});
                double *tempD = A->D[right-1];
                for(int row = 0;row < A->dSizes[right-1].first; row++)
                {
                    for(int col = 0; col < A->dSizes[right-1].second;col++)
                    {
                        tempD[col+(row*(A->dSizes[right-1].second))] -= (tempTTt[col+(row*(A->uSizes[right-1].first))] / B_c1_norm);
                    }
                }
                delete[] tempTTt;
                delete[] tempT;
                tempU = NULL;
                tempB = NULL;
                tempD = NULL;
            }

            else
            {
                double *tempU = A->U[right-1];
                double *tempUUt = new double[(A->uSizes[right-1].first)*(A->uSizes[right-1].first)];
                memset(tempUUt,0,sizeof(double)*(A->uSizes[right-1].first)*(A->uSizes[right-1].first));
             
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(A->uSizes[right-1].first),(A->uSizes[right-1].first),(A->uSizes[right-1].second),alpha,tempU,(A->uSizes[right-1].second),tempU,(A->uSizes[right-1].second),beta,tempUUt,(A->uSizes[right-1].first));
                double *tempD = A->D[right-1];
                for(int row = 0;row < A->dSizes[right-1].first; row++)
                {
                   for(int col = 0; col < A->dSizes[right-1].second;col++)
                   {
                        tempD[col+(row*(A->dSizes[right-1].second))] -= (tempUUt[col+(row*(A->uSizes[right-1].first))]*B_c1_norm);
                   }
                }                               
                delete [] tempUUt;
                tempU = NULL;
                tempD = NULL;
            }
        }
        else
        {
             //temporary pointer and pair for pushing onto the stack
             double *tempV;
             std::pair<int,int>tempP;
             if(A->bSizes[left-1].first <= A->bSizes[left-1].second)
             {
                double *Sp = new double[(A->bSizes[left-1].second)*(A->bSizes[left-1].first)];
                memset(Sp,0,sizeof(double)*(A->bSizes[left-1].second)*(A->bSizes[left-1].first));
                double *tempBC1 = new double[A->bSizes[left-1].first * (A->bSizes[left-1].second)];
                memcpy(tempBC1, A->B[left-1], sizeof(double)*(A->bSizes[left-1].first * (A->bSizes[left-1].second)));
                GetTransposeInPlace(tempBC1, A->bSizes[left-1].first, A->bSizes[left-1].second);
                //check if it can be made more efficient check weather it is corrrect or not?
                for(int row =0; row<A->bSizes[left-1].second;row++)
                {
                    for(int col=0;col<A->bSizes[left-1].first;col++)
                    {
                        Sp[col + row*(A->bSizes[left-1].first)] = tempBC1[col + row*(A->bSizes[left-1].first)] / sqrt_B_c1_norm;
                       // Sp[col+(row*(A->bSizes[left-1].first))] = (A->B[left-1][row+(col*(A->bSizes[left-1].first))])/sqrt_B_c1_norm;
                    }
                }
                delete[] tempBC1;
                tempV = Sp;                
                tempP = {(A->bSizes[left-1].second),(A->bSizes[left-1].first)};         
             }
             else
             {
                std::vector<int> child = bt->GetChildren(right);
                int eye_size = A->rSizes[(child[0]-1)].second;
                double *Sp   = new double[eye_size*eye_size];
                memset(Sp,0,sizeof(double)*eye_size*eye_size);
               
                for(int row_col = 0;row_col < eye_size; row_col++){
                    Sp[row_col+row_col*eye_size] = sqrt_B_c1_norm;
                }
                
                tempV = Sp;
                tempP = {eye_size,eye_size};
                
             }
//Note: check all the updates of the generator B
            std:: stack<pair<double*,pair<int,int>>> S;
            S.push({tempV,tempP});
            for(int j=right-1;j >= desc[right];j--)
            {               
                double *Sp = S.top().first;
                int Sp_row = (S.top().second).first;
                int Sp_col = (S.top().second).second;    
                S.pop();
                std::vector<int> temp_ch = bt->GetChildren(bt->tr[j-1]);
                
                if(j== temp_ch[1])
                {
                    int sib = temp_ch[0];
                   
                    S.push({Sp,{Sp_row,Sp_col}});

                    //Sp*Sp'
                    double *tempSST = new double[Sp_row*Sp_row];
                    memset(tempSST,0,sizeof(double)*Sp_row*Sp_row);                    
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,Sp_row,Sp_row,Sp_col,alpha,Sp,Sp_col,Sp,Sp_col,beta,tempSST,Sp_row);

                    //R{sib}*(Sp*Sp')
                    double *R_sib = A->R[(sib-1)]; 
                    double *tempRSST = new double[(A->rSizes[(sib-1)].first)*Sp_row];
                    memset(tempRSST,0,sizeof(double)*(A->rSizes[(sib-1)].first)*Sp_row);
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(A->rSizes[(sib-1)].first),Sp_row,Sp_row,alpha,R_sib,(A->rSizes[(sib-1)].second),tempSST,Sp_row,beta,tempRSST,Sp_row);

                    double *tempRSST_Rjt = new double[(A->rSizes[(sib-1)].first)*(A->rSizes[(j-1)].first)];
                    memset(tempRSST_Rjt,0,sizeof(double)*(A->rSizes[(sib-1)].first)*(A->rSizes[(j-1)].first));
                    double *Rj = A->R[j-1];
                   
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(A->rSizes[(sib-1)].first),(A->rSizes[(j-1)].first),(A->rSizes[(j-1)].second),alpha,tempRSST,Sp_row,Rj,(A->rSizes[(j-1)].second),beta,tempRSST_Rjt,(A->rSizes[(j-1)].first));

                    for(int row = 0; row<A->bSizes[sib-1].first;row++)
                    {
                        for(int col=0; col<A->bSizes[sib-1].second;col++)
                        {
                            A->B[sib-1][col+(row*(A->bSizes[sib-1].second))]=A->B[sib-1][col+(row*(A->bSizes[sib-1].second))]-tempRSST_Rjt[col+row*A->rSizes[j-1].first];
                        }
                    }

                    delete [] tempSST;
                    delete [] tempRSST;
                    delete [] tempRSST_Rjt;
                }
                //Sj = R{j} * Sp;
                double *Sj = new double[A->rSizes[j-1].first*Sp_col];
                //check this as well
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,A->rSizes[j-1].first,Sp_col,Sp_row,alpha,A->R[j-1],A->rSizes[j-1].second,Sp,Sp_col,beta,Sj,Sp_col); 
                std::vector<int> temp_ch2 = bt->GetChildren(j);
                if(temp_ch2.size()==0)
                {
                    //T = U{j}*Sj
                    double *tempT = new double[A->uSizes[j-1].first*Sp_col];
                    memset(tempT,0,sizeof(double)*(A->uSizes[j-1].first)*Sp_col);
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,A->uSizes[j-1].first,Sp_col,A->uSizes[j-1].second,alpha,A->U[j-1],A->uSizes[j-1].second,Sj,Sp_col,beta,tempT,Sp_col);                 
                    //T*T'
                    double *tempTTt = new double[(A->uSizes[j-1].first)*(A->uSizes[j-1].first)]; 
                    memset(tempTTt,0,sizeof(double)*(A->uSizes[j-1].first)*(A->uSizes[j-1].first));  
                    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,A->uSizes[j-1].first,A->uSizes[j-1].first,Sp_col,alpha,tempT,Sp_col,tempT,Sp_col,beta,tempTTt,A->uSizes[j-1].first);
                    //D{j} = D{j}-T*T' 
                    for(int row =0;row < (A->dSizes[j-1].first);row++)
                    {
                        for(int col=0;col < (A->dSizes[j-1].second);col++)
                        {
                            A->D[j-1][col+(row*(A->dSizes[j-1].second))] -= tempTTt[col+(row*(A->uSizes[j-1].first))];
                        }
                    }
                    delete [] tempT;
                    delete [] tempTTt;
                    delete [] Sj;
                }
                else
                {
                    std::vector<double>temp(Sj,Sj+(A->rSizes[j-1].first*Sp_col));
                    S.push({Sj,{A->rSizes[j-1].first,Sp_col}});
                }

              //  delete[] Sj;
              //  delete[] Sp;
            }
            //If stack is not empty
            while (!S.empty())
            {
                double *temp = S.top().first;
                S.pop();
                delete [] temp;
                temp=NULL;
            }       
        }

    }

    //Rank r-updates at non-leaf nodes
     //To store the diagonal blocks
    double** Z = new double*[n];

    //create a list to hold the matrix dimensions of the rank r updates 
	std::pair<int, int>* zSizes = new std::pair<int, int>[n];
    for(int i = 0; i < n; i++)
		zSizes[i]=std::make_pair(0,0); 
    
    //First double* is for array stored and int is for dimension 
    std:: stack<pair<double*,pair<int,int>>> SU;

    for(int i=0 ; i<n; i++)
    {
        std::vector<int> ch = bt->GetChildren(i+1);
        if(ch.empty())
        {
            SU.push(make_pair(A->U[i],A->uSizes[i]));            
            continue;
        }

        int left  = ch[0];
        int right = ch[1];        
        left = left - 1;
        right = right - 1;

        double *Uc2   = SU.top().first;
        std::pair<int,int>index_Uc2 = SU.top().second;
        SU.pop();
        
        double *Uc1   = SU.top().first;
        std::pair<int,int>index_Uc1 = SU.top().second;   
        SU.pop();
        
        double sqrt_B_c1_norm = 0;
        double B_c1_norm = 0;
        
        double *Temp_BC1 = new double[(A->bSizes[left].first)*(A->bSizes[left].second)];
        memcpy(Temp_BC1,A->B[left],sizeof(double) * (A->bSizes[left].first)*(A->bSizes[left].second));
        
        norm_svd(Temp_BC1, A->bSizes[left], &B_c1_norm);
        sqrt_B_c1_norm=std::sqrt(B_c1_norm);
        delete [] Temp_BC1;

        if(A->bSizes[left].first <= A->bSizes[left].second)
        {
            Z[i] = new double[(index_Uc1.first + index_Uc2.first)*(A->bSizes[left].first)];
            zSizes[i] = {(index_Uc1.first + index_Uc2.first),(A->bSizes[left].first)};
            memset(Z[i],0,sizeof(double)*(index_Uc1.first + index_Uc2.first)*(A->bSizes[left].first));
           
            memcpy(Z[i], Uc1, sizeof(double)*(index_Uc1.first*index_Uc1.second));
            for(int row_col = 0; row_col < ((index_Uc1.first*index_Uc1.second)); row_col++)
                Z[i][row_col] = Z[i][row_col]*sqrt_B_c1_norm;
           
            //Uc2*B{c1}'
            double *tempBC1 = new double[A->bSizes[left].first * (A->bSizes[left].second)];
            memcpy(tempBC1, A->B[left], sizeof(double)*(A->bSizes[left].first * (A->bSizes[left].second)));
            GetTransposeInPlace(tempBC1, A->bSizes[left].first, A->bSizes[left].second);

            double *temp_Uc2_Bc1 = new double[index_Uc2.first*(A->bSizes[left].first)];
            memset(temp_Uc2_Bc1, 0, sizeof(double)*index_Uc2.first*(A->bSizes[left].first));

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, index_Uc2.first, (A->bSizes[left].first), index_Uc2.second, alpha, Uc2, index_Uc2.second, tempBC1, A->bSizes[left].first, beta, temp_Uc2_Bc1, A->bSizes[left].first);

            for(int row_col = 0; row_col < (index_Uc2.first*(A->bSizes[left].first)); row_col++)
                temp_Uc2_Bc1[row_col] = temp_Uc2_Bc1[row_col] / sqrt_B_c1_norm;

            memcpy(Z[i]+(index_Uc1.first*index_Uc1.second), temp_Uc2_Bc1, sizeof(double) * (index_Uc2.first*(A->bSizes[left].first)));
            
            delete[] tempBC1;
        
            delete[] temp_Uc2_Bc1;
        }

        else
        {
            Z[i] = new double[(index_Uc1.first+index_Uc2.first)*(A->bSizes[left].second)];
            zSizes[i] = {(index_Uc1.first+index_Uc2.first),(A->bSizes[left].second)};
            memset(Z[i], 0, sizeof(double)*(index_Uc1.first+index_Uc2.first)*(A->bSizes[left].second));

            double *temp_Uc1_Bc1 = new double[index_Uc1.first*(A->bSizes[left].second)];
            memset(temp_Uc1_Bc1, 0, sizeof(double)*index_Uc1.first*(A->bSizes[left].second));
            
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,index_Uc1.first,(A->bSizes[left].second),index_Uc1.second,alpha,Uc1,index_Uc1.second,A->B[left],A->bSizes[left].second,beta,temp_Uc1_Bc1,A->bSizes[left].second);
            
            for(int row_col = 0; row_col < (index_Uc1.first * (A->bSizes[left].second)); row_col++)
                temp_Uc1_Bc1[row_col] = temp_Uc1_Bc1[row_col] / sqrt_B_c1_norm;
            
            memcpy(Z[i], temp_Uc1_Bc1, sizeof(double)*(index_Uc1.first * (A->bSizes[left].second)));
                      
            for(int row = index_Uc1.first, ind = 0; row < (index_Uc1.first+index_Uc2.first), ind < index_Uc2.first; row++, ind++)
            {
                for(int col = 0; col < index_Uc2.second; col++)
                {
                    Z[i][col + row * index_Uc2.second] = Uc2[col + ind * index_Uc2.second] * sqrt_B_c1_norm;
                }
            }
            delete [] temp_Uc1_Bc1;      
        }

        if(i < (n-1))
        {
            double *Ui = new double[(index_Uc1.first + index_Uc2.first) * (A->rSizes[left].second)];
            memset(Ui,0,sizeof(double)*(index_Uc1.first + index_Uc2.first) * (A->rSizes[left].second));
            
            double *Uc1_Rc1 = new double[index_Uc1.first * (A->rSizes[left].second)];
            memset(Uc1_Rc1, 0, sizeof(double)*index_Uc1.first * (A->rSizes[left].second));
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, index_Uc1.first, (A->rSizes[left].second), index_Uc1.second, alpha, Uc1, index_Uc1.second, A->R[left], A->rSizes[left].second, beta, Uc1_Rc1, (A->rSizes[left].second));

            double *Uc2_Rc2 = new double[index_Uc2.first * (A->rSizes[right].second)];
            memset(Uc2_Rc2,0,sizeof(double)*index_Uc2.first * (A->rSizes[right].second));
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,index_Uc2.first,(A->rSizes[right].second),index_Uc2.second,alpha,Uc2,index_Uc2.second,A->R[right],A->rSizes[right].second,beta,Uc2_Rc2,A->rSizes[right].second);

            // Ui = [Uc1 * R{c1}; Uc2 * R{c2}];
            memcpy(Ui,Uc1_Rc1,sizeof(double) * (index_Uc1.first * (A->rSizes[left].second)));
            memcpy(Ui+(index_Uc1.first * (A->rSizes[left].second)), Uc2_Rc2, sizeof(double) * index_Uc2.first*(A->rSizes[right].second));
            
            SU.push(make_pair(Ui,make_pair((index_Uc1.first + index_Uc2.first),(A->rSizes[left].second))));
            delete [] Uc2_Rc2;
            delete [] Uc1_Rc1;
        }

        delete [] Uc1;
        delete [] Uc2;
    }
    
    desc.clear();
    desc.shrink_to_fit();

    delete [] A->B;
    delete [] A->bSizes;
    delete [] A->R;
    delete [] A->rSizes;
    delete [] A->U;
    delete [] A->uSizes;
    ret->Z = Z;
    ret->zSizes = zSizes;
    ret->D = A->D;
    ret->dSizes = A->dSizes;
	return ret;
}