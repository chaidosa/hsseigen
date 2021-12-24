#include<string.h>
#include "BinTree.h"
#include "test_mat2hsssym.h"
#include "superDC.h"
#include "divide2.h"
#include "superdcmv_desc.h"
//#include "secualr"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}
using namespace std;

SDC* superDC(tHSSMat* A,  BinTree* bt, int* m, int mSize)
{

    //Dividing Stage
    DVD *resDvd = divide2(A,bt,m,mSize);

    //Conquering stage
    int N  = bt->GetNumNodes();
    //Index range of each node
    std::pair<int,int>* l = new std::pair<int,int>[N];

    for(int k = 0; k < N; k++)
		l[k] = std::make_pair(0,0);    

    l[0]                  = {0,m[0]-1};
    int it                = 0;
    int lt                = 0;

    for(int k = 0; k < N; k++)
    {
        std::vector<int> ch = bt->GetChildren(k+1);
        if(ch.size() == 0)
        {
            l[k] = {lt,lt+m[it]-1};
            lt   = l[k].second +1;
            it   = it+1;
        }
        else
        {
            l[k] = {l[ch[0]-1].first,l[ch[1]-1].second};

        }
    }
    /*
        Q0 = cell(k,1);         
        Lam = cell(k,1);
        rho = cell(k,1);
    */
    double **Q0 = new double*[N];    
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
        std::vector<int> ch = bt->GetChildren(i+1);
        //if ch is a leaf node
        if(ch.size() == 0){           
            //on exit of LAPACKE_dsyevd, it will have eigenvalues          
            double *d  = new double[resDvd->dSizes[i].first];
            double *e  = new double[resDvd->dSizes[i].first-1];
            double *t   = new double[resDvd->dSizes[i].first-1];
            //on exit of LAPACKE_dsyevd, it contains eigenvectors
            double *D   = new double[(resDvd->dSizes[i].first*(resDvd->dSizes[i].second))];
            memcpy(D,resDvd->D[i],sizeof(double)*(resDvd->dSizes[i].first)*(resDvd->dSizes[i].second));

            //computes all eigenvalues and eigenvectors of a real symmetric matrix D using divide and conquer algorithm       
            //by first converting matrix to tridiagonal form  
            int info = LAPACKE_dsytrd(LAPACK_ROW_MAJOR,'U',(resDvd->dSizes[i].first),D,(resDvd->dSizes[i].first),d,e,t);
            if(info > 0){
                cout<<"Error at converting";
                exit(1);
            }

            info = LAPACKE_dorgtr(LAPACK_ROW_MAJOR,'U',(resDvd->dSizes[i].first),D,(resDvd->dSizes[i].first),t);
            if(info > 0){
                cout<<"Error at converting";
                exit(1);
            } 

            info=LAPACKE_dstedc(LAPACK_ROW_MAJOR,'V',(resDvd->dSizes[i].first),d,e,D,(resDvd->dSizes[i].first));              
            //check for convergence
            if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
            }            
            
            Lam[i]      = d;
            LamSizes[i] = (resDvd->dSizes[i].first);
            Q0[i]       = D;            
            q0Sizes[i]  ={(resDvd->dSizes[i].first),(resDvd->dSizes[i].second)};  
            d = NULL;
            D = NULL;
            delete [] e;
            delete [] t;
        }        
        else
        {
            int left  = ch[0];
            int right = ch[1];

            superdcmv_desc(Q0,q0Sizes,resDvd->Z[i],resDvd->zSizes[i],bt,i,1,l,1024);
            //[Z{i}, nflops1] =  superdcmv_desc(Q0, Z{i}, tr, i, 1, rg, desc, N);
            /*
            Lam{i} = [Lam{c1}; Lam{c2}];
            Lam{c1} = [];
            Lam{c2} = [];
            r = size(Z{i}, 2);
            rho{i} = zeros(r, 1);
            for j = 1:r    
                [Lam{i}, Q0{i}{j}, nflops1, rho1] = secular(Lam{i}, Z{i}(:, j), tol, N);
                flops_conquer = flops_conquer + nflops1;
                rho{i}(j) = rho1;
                if j < r
                [Z{i}(:, j+1:r), nflops1] = superdcmv_cauchy(Q0{i}{j}, Z{i}(:, j+1:r), 1, N);  
                flops_conquer = flops_conquer + nflops1;
                end
            end
            */
            Lam[i] = new double[(LamSizes[left-1])+(LamSizes[right-1])];
            memcpy(Lam[i],Lam[left-1],sizeof(double)*(LamSizes[left-1]));
            memcpy(Lam[i]+(LamSizes[left-1]),Lam[right-1],sizeof(double)*(LamSizes[right-1]));
            delete [] Lam[left-1];
            delete [] Lam[right-1];
            LamSizes[left-1]  = 0;
            LamSizes[right-1] = 0;
            int r             = resDvd->zSizes[i].second;

            for(int j = 0; j < r; j++)
            {
             //secular(Lam{i}, Z{i}(:, j), tol, N)   
            } 

        }

       

    }
     
	return NULL;

} 
