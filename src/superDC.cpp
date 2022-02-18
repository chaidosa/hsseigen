#include<bits/stdc++.h>
#include<string.h>
#include "BinTree.h"
#include "test_mat2hsssym.h"
#include "superDC.h"
#include "divide2.h"
#include "superdcmv_desc.h"
#include "superdcmv_cauchy.h"
#include "eigenmatrix.h"
#include "secular.h"
#include "band2hss.h"
//#include "secualr"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}
using namespace std;

SDC* superDC(tHSSMat *A,  BinTree* bt, int* m, int mSize)
{

    cout<<"Reached superDC\n";
    //Dividing Stage
    DVD *resDvd = divide2(A,bt,m,mSize);

    cout<<"Sucess Divide\n";
    //Conquering stage
    int N  = bt->GetNumNodes();
    //Index range of each node
    std::pair<int,int>* l = new std::pair<int,int>[N];

    for(int k = 0; k < N; k++)
		l[k] = std::make_pair(0,0);    

    l[0]                  = {0,m[0] - 1};
    int it                = 0;
    int lt                = 0;

    for(int k = 0; k < N; k++)
    {
        std::vector<int> ch = bt->GetChildren(k + 1);
        if(ch.size() == 0)
        {
            l[k] = {lt,lt + m[it] - 1};
            lt   = l[k].second + 1;
            it   = it + 1;
        }
        else
        {
            l[k] = {l[ch[0] - 1].first,l[ch[1] - 1].second};

        }
    }
    /*
        Q0 = cell(k,1);         
        Lam = cell(k,1);
        rho = cell(k,1);
    */
    EIG_MAT **Q0 = new EIG_MAT*[N]; //Nikhil: [] not necessary (multiple places).
    //double **Q0 = new double*[N];    
    std::pair<int, int>* q0Sizes = new std::pair<int, int>[N]; 
    for(int k = 0; k < N; k++)
		q0Sizes[k]=std::make_pair(0,0);    
    
    double **Lam  = new double*[N]; 
    int *LamSizes = new int[N];
    

    //double **rho = new double*[N]; 
    //std::pair<int, int>* rhoSizes = new std::pair<int, int>[N];
    //for(int k = 0; k < N; k++)
	//	rhoSizes[k]=std::make_pair(0,0);


    for(int i = 0; i < N; i++)
    {
        std::vector<int> ch = bt->GetChildren(i + 1);
        //if current index i is a leaf node
        if(ch.size() == 0){
            double *E = new double[resDvd->dSizes[i].first];
            double *EV = new double[resDvd->dSizes[i].first*resDvd->dSizes[i].second];
            int *Isuppz = new int[2*resDvd->dSizes[i].second];
           // int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', resDvd->dSizes[i].first, resDvd->D[i], resDvd->dSizes[i].second, E);
            double abstol = 1.234e-27;
            int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', resDvd->dSizes[i].first, resDvd->D[i], resDvd->dSizes[i].second, NULL, NULL, NULL, NULL,abstol, &resDvd->dSizes[i].first, E, EV, resDvd->dSizes[i].second, Isuppz);
            if(info > 0){
                cout<<"Eigensolver doesn't work";
                exit(1);
            }

            Lam[i] = E;
            E = NULL;
            
            LamSizes[i] = resDvd->dSizes[i].first;
            
            Q0[i] = new EIG_MAT();
            Q0[i]->Q0_leaf = EV;
            EV = NULL;
            Q0[i]->Q0_nonleaf = NULL;
            q0Sizes[i] = {resDvd->dSizes[i].first, resDvd->dSizes[i].second};
            delete [] Isuppz;

        } 
        //current index is non-leaf node       
        else
        {
            int left  = ch[0];
            int right = ch[1];

            superdcmv_desc(Q0,q0Sizes,&(resDvd->Z[i]),resDvd->zSizes[i],bt,i,1,l,1024);
          
            Lam[i] = new double[(LamSizes[left-1]) + (LamSizes[right-1])];

            memcpy(Lam[i], Lam[left-1], sizeof(double) * (LamSizes[left-1]));
            memcpy(Lam[i]+(LamSizes[left - 1]), Lam[right - 1], sizeof(double) * (LamSizes[right - 1]));

            LamSizes[i] = (LamSizes[left - 1]) + (LamSizes[right - 1]);

            delete [] Lam[left - 1];
            delete [] Lam[right - 1];

            LamSizes[left - 1]  = 0;
            LamSizes[right - 1] = 0;

            int r             = resDvd->zSizes[i].second;

            Q0[i]             = new EIG_MAT();
            Q0[i]->Q0_leaf    = NULL;

            nonleaf **n_leaf    = new nonleaf*[r];
            

            double *temp_d   = Lam[i];
            int temp_d_size = LamSizes[i];

            for(int j = 0; j < r; j++){                
                //Z{:, j}                
                double *tempZ = new double[resDvd->zSizes[i].first];

                for(int row = 0; row < resDvd->zSizes[i].first; row++){
                    tempZ[row] = resDvd->Z[i][j + row * r];
                }

                SECU *res_sec = new SECU();
                res_sec = secular(temp_d, temp_d_size, tempZ, resDvd->zSizes[i].first, 1024);
                n_leaf[j] = res_sec->Q;


                delete [] temp_d;
                temp_d = res_sec->Lam;
                
                
                if(j < (r-1))
                {
                    double *tempZi = new double[resDvd->zSizes[i].first * (r - (j + 1))];

                    for(int row = 0; row < resDvd->zSizes[i].first; row++)
                        memcpy(tempZi + row*(r-(j+1)), resDvd->Z[i] + j + 1 + row * (resDvd->zSizes[i].second), sizeof(double) * (r - (j + 1)) );
                    
                    superdcmv_cauchy(&(n_leaf[j]), {1, 7}, &tempZi, {resDvd->zSizes[i].first, (r- ( j + 1)) }, 1);

                    for(int row = 0; row < resDvd->zSizes[i].first; row++)
                        memcpy(resDvd->Z[i] + j + 1 + row * (resDvd->zSizes[i].second), tempZi + row*(r - (j + 1)), sizeof(double) * (r - (j + 1)));
                    
                    delete [] tempZi;
                }
                
                delete [] tempZ;
                //will add rho later
            }
            Lam[i] = temp_d;
            Q0[i]->Q0_nonleaf = n_leaf;
            Q0[i]->n_non_leaf = r;
            q0Sizes[i] = {1, r};
        }

       

    }
    
    vector<double> tempeig;
    for(int k = 0; k < 256; k++)
        tempeig.push_back(Lam[N-1][k]);
    sort(tempeig.begin(), tempeig.end());
    for(int k = 0; k < 256; k++)
        cout<<setprecision(16)<<tempeig[k]<<endl;
	return NULL;

} 
