

#include<bits/stdc++.h>
#include<string.h>
#include "BinTree.h"
#include "superDC.h"
#include "divide2.h"
#include "superdcmv_desc.h"
#include "superdcmv_cauchy.h"
#include "eigenmatrix.h"
#include "secular.h"
#include "Generators.h"
#include "omp.h"
#include <sys/time.h>


#ifdef DIST
// extern "C"{
#include <mpi.h>
// }
#endif

extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}

//#include "/opt/intel/oneapi/mkl/2022.0.2/include/mkl_lapacke.h"
using namespace std;
std::pair<double *, nonleaf**> r_RankOneUpdate(double* Lam, int lamSize, std::pair<int, int>zSize, double* Z, nonleaf **n_leaf, int r);
std::pair<double*, double*> computeLeafEig(std::pair<int, int> dSize, double *D, int i);
void Eig_func(int i);

//Global declaration
int N;
DVD *resDvd;
EIG_MAT **Q0;
std::pair<int, int>* q0Sizes;
double **Lam; 
int *LamSizes;
BinTree *bt;
std::pair<int,int>* l;

void Eig_func(int i){
    vector<int> ch = bt->GetChildren(i + 1);
    
    if(ch.empty())
    {
        std::pair<double*, double *> E = computeLeafEig(make_pair(resDvd->dSizes[i].first, resDvd->dSizes[i].second), resDvd->D[i], i);
        Lam[i] = E.first;         
        LamSizes[i] = resDvd->dSizes[i].first;

        Q0[i] = new EIG_MAT();
        Q0[i]->Q0_leaf = E.second;            
        Q0[i]->Q0_nonleaf = NULL;
        q0Sizes[i] = {resDvd->dSizes[i].first, resDvd->dSizes[i].second};
       // light_switch[i] = true;
       return;
    }
    
    int left = ch[0];
    int right = ch[1];

    Eig_func(left-1);
    Eig_func(right-1);   
    
    //compute
    superdcmv_desc(Q0,q0Sizes,(resDvd->Z[i]),resDvd->zSizes[i],bt,i,1,l,fmmTrigger);           
    Lam[i] = new double[(LamSizes[left-1]) + (LamSizes[right-1])];

    std::copy(Lam[left-1], Lam[left-1] + LamSizes[left-1], Lam[i]);
    std::copy(Lam[right-1], Lam[right-1] + LamSizes[right-1], Lam[i] + LamSizes[left - 1]);
            
    LamSizes[i] = (LamSizes[left - 1]) + (LamSizes[right - 1]);
    //std::sort(Lam[i], Lam[i]+LamSizes[i]);

    delete [] Lam[left - 1];
    delete [] Lam[right - 1];

    LamSizes[left - 1]  = 0;
    LamSizes[right - 1] = 0;

    int r             = resDvd->zSizes[i].second;

    Q0[i]             = new EIG_MAT();
    Q0[i]->Q0_leaf    = NULL;

    nonleaf **n_leaf    = new nonleaf*[r];
            
    std::pair<double *, nonleaf**> result = r_RankOneUpdate(Lam[i], LamSizes[i], resDvd->zSizes[i], resDvd->Z[i], n_leaf, r);
    Lam[i] = result.first;
    Q0[i]->Q0_nonleaf = result.second;
    Q0[i]->n_non_leaf = r;
    q0Sizes[i] = {1, r};
    return;
}






std::pair<double *, nonleaf**> r_RankOneUpdate(double* Lam, int lamSize, std::pair<int, int>zSize, double* Z, nonleaf **n_leaf, int r){
    
    std::pair<double *, nonleaf**> result;
    
    double *temp_d   = Lam;
    int temp_d_size = lamSize;
    for(int j = 0; j < r; j++)
    {
                          
        double *tempZ = new double[zSize.first];

        for(int row = 0; row < zSize.first; row++){
            tempZ[row] = Z[j + row * r];
        }

        SECU *res_sec;
        res_sec = secular(temp_d, temp_d_size, tempZ, zSize.first,  fmmTrigger);
        n_leaf[j] = res_sec->Q;
                
        delete [] temp_d;
        temp_d = res_sec->Lam;
                
                
        if(j < (r-1))
        {
            double *tempZi = new double[zSize.first * (r - (j + 1))];

            for(int row = 0; row < zSize.first; row++)
                memcpy(tempZi + row*(r-(j+1)), Z + (j + 1) + row * (zSize.second), sizeof(double) * (r - (j + 1)) );
                    
            tempZi = superdcmv_cauchy((n_leaf[j]), {1, 7}, tempZi, {zSize.first, (r- ( j + 1)) }, 1, fmmTrigger);

            for(int row = 0; row < zSize.first; row++)
                memcpy(Z + j + 1 + row * (zSize.second), tempZi + row*(r - (j + 1)), sizeof(double) * (r - (j + 1)));
                    
            delete [] tempZi;
        }
                
    delete [] tempZ;
                //will add rho later
    }
    
    result.first = temp_d;
    result.second = n_leaf;
    return result;       
}



std::pair<double*, double*> computeLeafEig(std::pair<int, int> dSize, double *D, int i){
            
    std::pair<double*, double*> result;
    double *E = new double[dSize.first];
    double *EV = new double[dSize.first*dSize.second];
    int *Isuppz = new int[2*dSize.second];
            
    double abstol = 1.234e-27;
        
    //cout << "Computing eigenvalues and eigenvectors for node: "<<(i+1)<<"\n";
    int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', dSize.first, D, dSize.second, 0., 0., 0, 0,abstol, &dSize.first, E, EV, dSize.second, Isuppz);
            
    if(info > 0){
        cout<<"Eigensolver doesn't work";
        exit(1);
    }
            
    result.first = E;
    result.second = EV;
    E = NULL;        
           
    delete [] Isuppz;
    return result;
}


SDC* superDC(GEN *A, BinTree* btree, int* m, int mSize, int nProc)
{

    struct timeval timeStart, timeEnd;
    gettimeofday(&timeStart, 0);

      
bt = btree;
   cout<<"Reached superDC\n";
    //Dividing Stage
    resDvd = divide2(A, bt, m, mSize);
    
   cout<<"Success Divide\n";
    //Conquering stage
    N  = bt->GetNumNodes();

    Q0 = new EIG_MAT*[N];
    q0Sizes = new std::pair<int, int>[N];
    Lam  = new double*[N]; 
    LamSizes = new int[N];
    l = new std::pair<int,int>[N];  

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
 
    for(int k = 0; k < N; k++)
		q0Sizes[k]=std::make_pair(0,0);    
 
   // vector<int> counter(N+1, 0);
   // std::vector<int> WorkQueue(bt->leaves.begin(), bt->leaves.end());
    

  //  struct timeval timeStart, timeEnd;
   // gettimeofday(&timeStart, 0);
   // cout<<"Number of processors:"<<omp_get_num_procs()<<endl;


#if defined(PARALLEL)
    vector<int> counter(N+1, 0);
    std::vector<int> WorkQueue(bt->leaves.begin(), bt->leaves.end());

   // cout<<"Number of processors:"<<omp_get_num_procs()<<endl;

    omp_set_num_threads(nProc);
    #pragma omp parallel
    {

        //omp_set_num_threads(omp_get_num_procs());
       // cout<<"In Parallel:"<<omp_in_parallel()<<endl;
      //  cout<<"Number of threads:"<<omp_get_num_threads()<<endl;
        #pragma omp for
        for(int thr = 0; thr < WorkQueue.size(); thr++){
          // cout<<"iteration:"<<thr<<"thread:"<<omp_get_thread_num()<<endl;            
            while(true){
                int node;
                vector<int> ch;
                #pragma omp critical
                {
                    if(!WorkQueue.empty()){
                        node = WorkQueue.back();
                        WorkQueue.pop_back();
                    }
                    else{
                        node = -1;
                    }
                }

                if(node != -1 && bt->GetChildren(node).empty()){
                    int i = node - 1;
                    std::pair<double*, double *> E = computeLeafEig(make_pair(resDvd->dSizes[i].first, resDvd->dSizes[i].second), resDvd->D[i], i);
                    Lam[i] = E.first;
                    LamSizes[i] = resDvd->dSizes[i].first;

                                    Q0[i] = new EIG_MAT();
                    Q0[i]->Q0_leaf = E.second;
                    Q0[i]->Q0_nonleaf = NULL;
                    q0Sizes[i] = {resDvd->dSizes[i].first, resDvd->dSizes[i].second};

                    #pragma omp critical
                    {
                        counter[bt->tr[node - 1]]++;
                        if(counter[bt->tr[node - 1]] == 2){
                                WorkQueue.push_back(bt->tr[node - 1]);
                            }
                    }
                }

                else if (node != -1 && !bt->GetChildren(node).empty()){
                   // cout<<"Computing Internal node: "<<node<<"\n";
                    ch = bt->GetChildren(node);
                    int left = ch[0];
                    int right = ch[1];
                    int i = node - 1;
                    resDvd->Z[i] = superdcmv_desc(Q0,q0Sizes,(resDvd->Z[i]),resDvd->zSizes[i],bt,i,1,l,fmmTrigger);
                    Lam[i] = new double[(LamSizes[left-1]) + (LamSizes[right-1])];

                    std::copy(Lam[left-1], Lam[left-1] + LamSizes[left-1], Lam[i]);
                    std::copy(Lam[right-1], Lam[right-1] + LamSizes[right-1], Lam[i] + LamSizes[left - 1]);

                    LamSizes[i] = (LamSizes[left - 1]) + (LamSizes[right - 1]);
                    //std::sort(Lam[i], Lam[i]+LamSizes[i]);

                    delete [] Lam[left - 1];
                    delete [] Lam[right - 1];

                    LamSizes[left - 1]  = 0;
                    LamSizes[right - 1] = 0;

                    int r             = resDvd->zSizes[i].second;

                    Q0[i]             = new EIG_MAT();
                    Q0[i]->Q0_leaf    = NULL;

                    nonleaf **n_leaf    = new nonleaf*[r];


                    std::pair<double *, nonleaf**> result = r_RankOneUpdate(Lam[i], LamSizes[i], resDvd->zSizes[i], resDvd->Z[i], n_leaf, r);
                    Lam[i] = result.first;
                    Q0[i]->Q0_nonleaf = result.second;
                    Q0[i]->n_non_leaf = r;
                    q0Sizes[i] = {1, r};

                    #pragma omp critical
                    {
                        counter[bt->tr[node - 1]]++;
                        if(counter[bt->tr[node - 1]] == 2){
                                WorkQueue.push_back(bt->tr[node - 1]);
                        }

                    }
                }
                else if(node == -1)
                    break;
            }

        }
    }    
#else
  std::queue<int> Work;
  std::vector<int> counter(N+1, 0);

  for(int itr = 0; itr < bt->leaves.size(); itr++){
        Work.push(bt->leaves[itr]);
  }

  while(!Work.empty()){
        vector<int> ch;
        int node = Work.front();
        if(node != N)
            counter[bt->tr[node-1]]++;
        Work.pop();
        if(counter[bt->tr[node-1]] == 2)
          Work.push(bt->tr[node-1]);

        if(bt->GetChildren(node).empty()){
                    int i = node - 1;
                    std::pair<double*, double *> E = computeLeafEig(make_pair(resDvd->dSizes[i].first, resDvd->dSizes[i].second), resDvd->D[i], i);

    
               Lam[i] = E.first;
                    LamSizes[i] = resDvd->dSizes[i].first;

                    Q0[i] = new EIG_MAT();
                    Q0[i]->Q0_leaf = E.second;
                    Q0[i]->Q0_nonleaf = NULL;
                    q0Sizes[i] = {resDvd->dSizes[i].first, resDvd->dSizes[i].second};

        }

        else if(!bt->GetChildren(node).empty()){
                ch = bt->GetChildren(node);
                int left = ch[0];
                int right = ch[1];
                int i = node - 1;
                resDvd->Z[i] = superdcmv_desc(Q0,q0Sizes,(resDvd->Z[i]),resDvd->zSizes[i],bt,i,1,l,fmmTrigger);
                Lam[i] = new double[(LamSizes[left-1]) + (LamSizes[right-1])];
                std::copy(Lam[left-1], Lam[left-1] + LamSizes[left-1], Lam[i]);
                std::copy(Lam[right-1], Lam[right-1] + LamSizes[right-1], Lam[i] + LamSizes[left - 1]);

                LamSizes[i] = (LamSizes[left - 1]) + (LamSizes[right - 1]);
                delete [] Lam[left - 1];
                delete [] Lam[right - 1];

                LamSizes[left - 1] = 0;
                LamSizes[right- 1] = 0;

                int r = resDvd->zSizes[i].second;
                    // cout << "R is "<<r<<"\n";

                Q0[i] = new EIG_MAT();
                Q0[i]->Q0_leaf = NULL;

                nonleaf **n_leaf = new nonleaf*[r];
                std::pair<double*, nonleaf**> result = r_RankOneUpdate(Lam[i], LamSizes[i], resDvd->zSizes[i], resDvd->Z[i], n_leaf, r);
                Lam[i] = result.first;
                Q0[i]->Q0_nonleaf = result.second;
                Q0[i]->n_non_leaf = r;
                q0Sizes[i] = {1, r};
        }
}
#endif

 
    gettimeofday(&timeEnd, 0);
    long long elapsed = (timeEnd.tv_sec-timeStart.tv_sec)*1000000LL + timeEnd.tv_usec-timeStart.tv_usec;
        printf ("\nDone. %f usecs\n",elapsed/(double)1000000);

    
    std::ofstream txtOut;
    txtOut.open("output.txt", std::ofstream::out);
    txtOut <<setprecision(10)<<elapsed/(double)1000000<<" seconds"<<endl;
    vector<double> tempeig;
    for(int k = 0; k < LamSizes[N-1]; k++)
        tempeig.push_back(Lam[N-1][k]);

    std::sort(tempeig.begin(), tempeig.end());
    int count = 0;
    for(int k = 0; k < LamSizes[N-1]; k++){
        count++;
        txtOut<<setprecision(20)<<tempeig[k]<<endl;
    }
    cout << count;
    SDC *resSDC = new SDC();

    resSDC->Q = Q0;
    resSDC->L = Lam[N-1];
    resSDC->lSize = LamSizes[N-1];
    resSDC->qSizes = q0Sizes;
    return resSDC;
}

#ifdef DIST
#include "distributedRoutines.h"


SDC *dsuperDc(GEN *A, BinTree *btree, int *m, int mSize, int nProc, MPI_Comm process_grid)
{
    // struct timeval timeStart, timeEnd;
    // gettimeofday(&timeStart, 0);
    bt = btree;
    // cout << "Reached superDC\n";
    // Dividing Stage
    resDvd = divide2(A, bt, m, mSize);

    // cout << "Success Divide\n";
    // Conquering stage
    N = bt->GetNumNodes();

    Q0 = new EIG_MAT *[N];
    q0Sizes = new std::pair<int, int>[N];
    Lam = new double *[N];
    LamSizes = new int[N];
    l = new std::pair<int, int>[N];

    for (int k = 0; k < N; k++)
        l[k] = std::make_pair(0, 0);

    l[0] = {0, m[0] - 1};
    int it = 0;
    int lt = 0;

    for (int k = 0; k < N; k++)
    {
        std::vector<int> ch = bt->GetChildren(k + 1);
        if (ch.size() == 0)
        {
            l[k] = {lt, lt + m[it] - 1};
            lt = l[k].second + 1;
            it = it + 1;
        }
        else
        {
            l[k] = {l[ch[0] - 1].first, l[ch[1] - 1].second};
        }
    }

    for (int k = 0; k < N; k++)
        q0Sizes[k] = std::make_pair(0, 0);

#if defined(DIST)

    /**
     * Get the grid information
    */
    int nprocs, myrank;
    MPI_Comm_size(process_grid, &nprocs);
    MPI_Comm_rank(process_grid, &myrank);

    

    // int coords[2];
    // int dims[2];
    // int reorder[2];
    // MPI_Cart_get(process_grid, 2, dims, reorder, coords);
    // int my_col = coords[1]; // for row wise ordering
    // int my_row = coords[0];

    // Timer start
    double tstart, tend;
    tstart = MPI_Wtime();

    /**
     * Compute the leaf eigenvalues, define communicators and distribute
     */

    vector<vector<int>> level_order_nodes = bt->nodeAtLvl;

    vector<vector<int>> core_map(level_order_nodes.size());
    vector<vector<int>> level_map(level_order_nodes.size());
    core_map[0].push_back(0);
    level_map[0].push_back(0);
    if (level_order_nodes.size() >= 2)
    {   
        core_map[1].resize(2);
        level_map[1].resize(2);
        core_map[1][0] = 0;
        core_map[1][1] = 1;

        level_map[1][0] = 0;
        level_map[1][1] = 1;
    }

    for (int level = 2; level < level_order_nodes.size(); level++)
    {
        int num_nodes = (int)pow(2, level);
        core_map[level].resize(num_nodes);
        level_map[level].resize(num_nodes);
        
        int even_step = (int)pow(2, level-1) -2;
        for(int i=0; i<num_nodes/2; i++){
            if (i%2==0)
            {
                core_map[level][i] = core_map[level-1][i/2];
                level_map[level][core_map[level-1][i/2]] = i; 
            } else
            {
                core_map[level][i] = even_step+2;
                level_map[level][even_step+2] = i;
                even_step+=2;
            }
        }

        int odd_step = (int)pow(2, level-1) -1;
        for (int i = num_nodes/2; i < num_nodes; i++)
        {
            if (i%2==0)
            {
                core_map[level][i] = core_map[level-1][i/2];
                level_map[level][core_map[level-1][i/2]] = i; 
            } else
            {
                core_map[level][i] = odd_step+2;
                level_map[level][odd_step+2] = i;
                odd_step+=2;
            }
        }
    }

    // print out level order tree and core map
    // printHssTree(level_order_nodes, core_map);

    // A set to store the Qo[indx] which the process has
    set<int> Q0_list;

    for (int level = level_order_nodes.size() - 1; level >= 0; level--)
    {   

        vector<int> work_done;
        for (int node_index = 0; node_index < level_order_nodes[level].size(); node_index++)
        {
            if (core_map[level][node_index] == myrank)
            {   
                work_done.push_back(node_index);
                int node = level_order_nodes[level][node_index];
                if (bt->GetChildren(node).empty())
                {
                    int i = node - 1;
                    std::pair<double *, double *> E = computeLeafEig(make_pair(resDvd->dSizes[i].first, resDvd->dSizes[i].second), resDvd->D[i], i);

                    Lam[i] = E.first;
                    LamSizes[i] = resDvd->dSizes[i].first;
                    
                    Q0[i] = new EIG_MAT();
                    Q0[i]->Q0_leaf = E.second;
                    Q0[i]->Q0_nonleaf = NULL;
                    q0Sizes[i] = {resDvd->dSizes[i].first, resDvd->dSizes[i].second};
                    
                    Q0_list.insert(i);
                }

                else if (!bt->GetChildren(node).empty())
                {
                    vector<int> ch = bt->GetChildren(node);
                    int left = ch[0];
                    int right = ch[1];
                    int i = node - 1;
                    resDvd->Z[i] = superdcmv_desc(Q0, q0Sizes, (resDvd->Z[i]), resDvd->zSizes[i], bt, i, 1, l, fmmTrigger);

                    Lam[i] = new double[(LamSizes[left - 1]) + (LamSizes[right - 1])];
                    std::copy(Lam[left - 1], Lam[left - 1] + LamSizes[left - 1], Lam[i]);
                    std::copy(Lam[right - 1], Lam[right - 1] + LamSizes[right - 1], Lam[i] + LamSizes[left - 1]);

                    LamSizes[i] = (LamSizes[left - 1]) + (LamSizes[right - 1]);
                    delete[] Lam[left - 1];
                    delete[] Lam[right - 1];

                    LamSizes[left - 1] = 0;
                    LamSizes[right - 1] = 0;

                    int r = resDvd->zSizes[i].second;

                    Q0[i] = new EIG_MAT();
                    Q0[i]->Q0_leaf = NULL;

                    nonleaf **n_leaf = new nonleaf *[r];
                    std::pair<double *, nonleaf **> result = r_RankOneUpdate(Lam[i], LamSizes[i], resDvd->zSizes[i], resDvd->Z[i], n_leaf, r);
                    Lam[i] = result.first;
                    Q0[i]->Q0_nonleaf = result.second;
                    Q0[i]->n_non_leaf = r;
                    q0Sizes[i] = {1, r};

                    Q0_list.insert(i);
                }
            }
        }

        int SEND_Error;
        if (myrank >= (int)pow(2, level - 1) && myrank < (int)pow(2, level) && level>0)
        { // send
            int myindx = level_map[level][myrank];
            int send_to = core_map[level][myindx - 1];

            // first send the number of results to sent
            int t = work_done.size();

            // first send contains sizes and if sending leaf or node
            SEND_Error= MPI_Send((void *)&t, 1, MPI_INT, send_to, 0, process_grid);
            if (SEND_Error!=MPI_SUCCESS)
            {
                cout << "SendError at 670" << "\n";
            }
            

            for (int i = 0; i < t; i++)
            {
                int node = level_order_nodes[level][work_done[i]];
                const int indx = node - 1;
                int LamSize = LamSizes[indx];
                int q0_size_first = q0Sizes[indx].first;
                int q0_size_second = q0Sizes[indx].second;
                bool is_leaf = bt->GetChildren(node).empty();
                double *send_buf = new double[6];

                send_buf[0] = indx;
                send_buf[1] = LamSize;
                send_buf[2] = q0_size_first;
                send_buf[3] = q0_size_second;
                send_buf[4] = is_leaf;
                send_buf[5] = Q0_list.size();

                SEND_Error= MPI_Send(send_buf, 6, MPI_DOUBLE, send_to, 0, process_grid);
                delete[] send_buf;                
                if (SEND_Error != MPI_SUCCESS)
                {
                    cout << "SendError at 696"
                         << "\n";
                }                        

                if (is_leaf)
                {

                    SEND_Error = MPI_Send(Lam[indx], LamSize, MPI_DOUBLE, send_to, 0, process_grid);
                    if (SEND_Error != MPI_SUCCESS)
                    {
                        cout << "SendError at 700"
                             << "\n";
                    }

                    SEND_Error = MPI_Send(Q0[indx]->Q0_leaf, q0_size_first * q0_size_second, MPI_DOUBLE, send_to, 0, process_grid);
                    if (SEND_Error != MPI_SUCCESS)
                    {
                        cout << "SendError at 670"
                             << "\n";
                    }
                }
                else
                {   
                    SEND_Error = MPI_Send(Lam[indx], LamSize, MPI_DOUBLE, send_to, 0, process_grid);
                    if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Lam at 700, rank: " << myrank << "\n";}

                    int* Q0_list_buff = new int[2*Q0_list.size()]; int step=0;
                    for (int x : Q0_list){
                        Q0_list_buff[step] = x;
                        step++;
                    }
                    for (int x : Q0_list){
                        Q0_list_buff[step] = Q0[x]->n_non_leaf;
                        step++;
                    }
                    SEND_Error = MPI_Send(Q0_list_buff, 2*Q0_list.size(), MPI_INT, send_to, 0, process_grid);
                    if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Q0_list at 700, rank: " << myrank << "\n";}
                    
                    

                    for (int Q0_indx : Q0_list){

                        int* Q0_indx_sizes = new int[4+Q0[Q0_indx]->n_non_leaf*(12+10+4)+2];
                        serializeQsizes(Q0[Q0_indx], Q0_indx_sizes, bt->GetChildren(Q0_indx+1).empty(), q0Sizes[Q0_indx]);
                        SEND_Error = MPI_Send(Q0_indx_sizes, 4+Q0[Q0_indx]->n_non_leaf*(12+10+4)+2, MPI_INT, send_to, 0, process_grid);
                        if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Q0_sizes at 700, rank: " << myrank << "\n";}

                        if (bt->GetChildren(Q0_indx+1).empty())
                        {
                            SEND_Error = MPI_Send(Q0[Q0_indx]->Q0_leaf, q0Sizes[Q0_indx].first*q0Sizes[Q0_indx].second, MPI_DOUBLE, send_to, 0, process_grid);
                            if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Q0_leaf at 700, rank: " << myrank << "\n";}
                        } 
                        else
                        {   
                            
                            int int_buff_sz = Q0_indx_sizes[4+Q0[Q0_indx]->n_non_leaf*(12+10+4)];
                            int double_buff_sz = Q0_indx_sizes[4+Q0[Q0_indx]->n_non_leaf*(12+10+4)+1];

                            int* T_ORG = new int[int_buff_sz];
                            double* double_buff = new double[double_buff_sz];

                            serializeQNonLeaf(Q0[Q0_indx], T_ORG, double_buff);

                            SEND_Error = MPI_Send(T_ORG, int_buff_sz, MPI_INT, send_to, 0, process_grid);
                            if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Lam at 700, rank: " << myrank << "\n";}

                            SEND_Error = MPI_Send(double_buff, double_buff_sz, MPI_DOUBLE, send_to, 0, process_grid);
                            if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Lam at 700, rank: " << myrank << "\n";}
                            // cout << "Sending sizes " << int_buff_sz << ", " << double_buff_sz << " myrank " << myrank << " sends to " << send_to << "\n";


                            // delete[] T_ORG;
                            // delete[] double_buff;
                        }

                        // delete[] Q0_indx_sizes;
                    }
                    
                    // delete[] Q0_list_buff;
                }
            }
        }
        else
        { // recv

            if (myrank < (int)pow(2, level) && level>0)
            {

                int myindx = level_map[level][myrank];
                int recv_from = core_map[level][myindx + 1];

                // number of results to recv
                int t = 0;
                MPI_Status status;
                MPI_Recv((void *)&t, 1, MPI_INT, recv_from, 0, process_grid, &status);

                for (int i = 0; i < t; i++)
                {
                    double *recv_sizes_first = new double[6];
                    MPI_Recv((void *)recv_sizes_first, 6, MPI_DOUBLE, recv_from, 0, process_grid, &status);

                    int indx = (int)recv_sizes_first[0];
                    int lamsz = (int)recv_sizes_first[1];
                    int qsz_first = (int)recv_sizes_first[2];
                    int qsz_second = (int)recv_sizes_first[3];
                    int isLeaf = (int)recv_sizes_first[4];
                    int Q0_list_size = (int)recv_sizes_first[5];


                    q0Sizes[indx] = make_pair(qsz_first, qsz_second);
                    LamSizes[indx] = lamsz;

                    // cout << "Myrank, recvd Qi: " << myrank << "," << indx << "\n";
                    if (isLeaf)
                    {
                        Lam[indx] = new double[lamsz];
                        Q0[indx] = new EIG_MAT();
                        Q0[indx]->Q0_nonleaf = NULL;
                        Q0[indx]->Q0_leaf = new double[qsz_first * qsz_first];
                        MPI_Recv(Lam[indx], lamsz, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                        MPI_Recv(Q0[indx]->Q0_leaf, qsz_first * qsz_first, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                        Q0_list.insert(indx);
                    }
                    else
                    {

                        Lam[indx] = new double[lamsz];
                        MPI_Recv(Lam[indx], lamsz, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                        

                        int* Q0_list_buff = new int[2*Q0_list_size];
                        MPI_Recv(Q0_list_buff, 2*Q0_list_size, MPI_INT, recv_from, 0, MPI_COMM_WORLD, &status);

                        

                        for (int _indx = 0; _indx < Q0_list_size; _indx++)
                        {   
                            int step = Q0_list_size;
                            int Q0_n_non_leaf = Q0_list_buff[step+_indx];
                            int Q0_indx = Q0_list_buff[_indx];
                            Q0[Q0_indx] = new EIG_MAT();
                            Q0[Q0_indx]->n_non_leaf = Q0_n_non_leaf;
                            int* Q0_sizes_buff = new int[4+Q0_n_non_leaf*(12+10+4)+2]; // error here
                            MPI_Recv(Q0_sizes_buff, 4+Q0_n_non_leaf*(12+10+4)+2, MPI_INT, recv_from, 0, process_grid, &status);
                            bool isleaf = Q0_sizes_buff[2];


                            

                            pair<int,int> buff_sizes = deserializeQsizes(Q0[Q0_indx], Q0_sizes_buff);
                            q0Sizes[Q0_indx] = make_pair(Q0_sizes_buff[0], Q0_sizes_buff[1]);
                            if (isleaf)
                            {
                                Q0[Q0_indx]->Q0_leaf = new double[q0Sizes[Q0_indx].first*q0Sizes[Q0_indx].second];
                                MPI_Recv(Q0[Q0_indx]->Q0_leaf, q0Sizes[Q0_indx].first*q0Sizes[Q0_indx].second, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                            }
                            else
                            {
                                int int_buff_sz = buff_sizes.first;
                                int double_buff_sz = buff_sizes.second;

                                int* T_ORG = new int[int_buff_sz];
                                double* double_buff = new double[double_buff_sz];

                                MPI_Recv(T_ORG, int_buff_sz, MPI_INT, recv_from, 0, process_grid, &status);

                                MPI_Recv(double_buff, double_buff_sz, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                                // if (level==1&&myrank==0)
                                // cout << "Reached Here 825 " << myrank << " " << Q0_indx << " " << Q0_sizes_buff[2] << " " << Q0_sizes_buff[3] << " " << qsz_second << "\n";
                                deserialzeQNonLeaf(Q0[Q0_indx], T_ORG, double_buff);
                                
                            }

                            Q0_list.insert(Q0_indx);
                            delete[] Q0_sizes_buff;
                        }
                        
                        delete[] Q0_list_buff;
                    }
                }
            }
        }

        // MPI_Barrier(process_grid);
        // cout << "Completed level " << level << " " << myrank << "\n";
        // if(myrank==0){
        //     cout << "\n";
        // }
        // sleep(1);
        MPI_Barrier(process_grid);
    }

    /**
     * Calculate execution time
     */
    tend = MPI_Wtime();
    double etime = tend - tstart;
    double max_etime = 0.0;
    MPI_Reduce(&etime, &max_etime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0)
    {
        cout << "Distributed Superdc took " << max_etime << " seconds"
             << "\n";
    }
    MPI_Barrier(process_grid);
    // MPI_Abort(MPI_COMM_WORLD, 0);
#endif

    // gettimeofday(&timeEnd, 0);
    // long long elapsed = (timeEnd.tv_sec - timeStart.tv_sec) * 1000000LL + timeEnd.tv_usec - timeStart.tv_usec;
    // printf("\nDone. %f usecs\n", elapsed / (double)1000000);
#ifdef DIST
    if (myrank==0){
#endif
    std::ofstream txtOut;
    txtOut.open("output.txt", std::ofstream::out);
    // txtOut << setprecision(10) << elapsed / (double)1000000 << " seconds" << endl;
    vector<double> tempeig;
    for (int k = 0; k < LamSizes[N - 1]; k++)
        tempeig.push_back(Lam[N - 1][k]);

    std::sort(tempeig.begin(), tempeig.end());
    int count = 0;
    for (int k = 0; k < LamSizes[N - 1]; k++)
    {
        count++;
        txtOut << setprecision(20) << tempeig[k] << endl;
    }
    cout << count;
    SDC *resSDC = new SDC();

    resSDC->Q = Q0;
    resSDC->L = Lam[N - 1];
    resSDC->lSize = LamSizes[N - 1];
    resSDC->qSizes = q0Sizes;
    return resSDC;
#ifdef DIST
    } else
    {
        return new SDC();
    }
    
#endif

}
#endif
/*
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
#include "Generators.h"
#include <sys/time.h>
//#include "secualr"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}
using namespace std;

SDC* superDC(GEN *A,  BinTree* bt, int* m, int mSize, int nProc)
{
    struct timeval timeStart, timeEnd;
    gettimeofday(&timeStart, 0);

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
    
        Q0 = cell(k,1);         
        Lam = cell(k,1);
        rho = cell(k,1);
    
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
        if(ch.size() == 0)
        {
            double *E = new double[resDvd->dSizes[i].first];
            double *EV = new double[resDvd->dSizes[i].first*resDvd->dSizes[i].second];
            int *Isuppz = new int[2*resDvd->dSizes[i].second];
           // int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', resDvd->dSizes[i].first, resDvd->D[i], resDvd->dSizes[i].second, E);
            double abstol = 1.234e-27;
        
            cout << "Computing eigenvalues and eigenvectors for node: "<<(i+1)<<"\n";
            int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', resDvd->dSizes[i].first, resDvd->D[i], resDvd->dSizes[i].second, NULL, NULL, NULL, NULL,abstol, &resDvd->dSizes[i].first, E, EV, resDvd->dSizes[i].second, Isuppz);
           // int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', resDvd->dSizes[i].first, resDvd->D[i], resDvd->dSizes[i].second, E);
             //int info = LAPACKE_dggev(LAPACK_ROW_MAJOR, 'V', 'N', resDvd->dSizes[i].first, resDvd->Z[i], resDvd->dSizes[i].second, )

            if(info > 0){
                cout<<"Eigensolver doesn't work";
                exit(1);
            }
            //EV = resDvd->D[i];
            std::sort(E, E+resDvd->dSizes[i].first);
            vector<double>ee(E, E+resDvd->dSizes[i].first);

            //std::copy(E, E+resDvd->dSizes[i].first, ee);
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
            //if(i == 2)
            //    continue;

            int left  = ch[0];
            int right = ch[1];

            superdcmv_desc(Q0,q0Sizes,&(resDvd->Z[i]),resDvd->zSizes[i],bt,i,1,l,1024);

            //vector<double>Zi(resDvd->Z[i], resDvd->Z[i]+(resDvd->zSizes[i].first * resDvd->zSizes[i].second));

            Lam[i] = new double[(LamSizes[left-1]) + (LamSizes[right-1])];

            std::copy(Lam[left-1], Lam[left-1] + LamSizes[left-1], Lam[i]);
            std::copy(Lam[right-1], Lam[right-1] + LamSizes[right-1], Lam[i] + LamSizes[left - 1]);
            
            LamSizes[i] = (LamSizes[left - 1]) + (LamSizes[right - 1]);
            //std::sort(Lam[i], Lam[i]+LamSizes[i]);

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

            for(int j = 0; j < r; j++)
            {                
                //Z{:, j}                
                double *tempZ = new double[resDvd->zSizes[i].first];

                for(int row = 0; row < resDvd->zSizes[i].first; row++)
                {
                    tempZ[row] = resDvd->Z[i][j + row * r];
                }

                SECU *res_sec;
                res_sec = secular(temp_d, temp_d_size, tempZ, resDvd->zSizes[i].first, 1024);
                n_leaf[j] = res_sec->Q;
                vector<double> LamVec(res_sec->Lam, (res_sec->Lam)+ temp_d_size);

                delete [] temp_d;
                temp_d = res_sec->Lam;
                
                
                if(j < (r-1))
                {
                    double *tempZi = new double[resDvd->zSizes[i].first * (r - (j + 1))];

                    for(int row = 0; row < resDvd->zSizes[i].first; row++)
                        memcpy(tempZi + row*(r-(j+1)), resDvd->Z[i] + (j + 1) + row * (resDvd->zSizes[i].second), sizeof(double) * (r - (j + 1)) );
                    
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

    gettimeofday(&timeEnd, 0);
    long long elapsed = (timeEnd.tv_sec-timeStart.tv_sec)*1000000LL + timeEnd.tv_usec-timeStart.tv_usec;
        printf ("\nDone. %f usecs\n",elapsed/(double)1000000);
    
        vector<double> tempeig;
        for(int k = 0; k < LamSizes[N-1]; k++)
            tempeig.push_back(Lam[N-1][k]);

        sort(tempeig.begin(), tempeig.end());
        int count = 0;
        for(int k = 0; k < LamSizes[N-1]; k++){
            count++;
            cout<<setprecision(16)<<tempeig[k]<<endl;
        }
	    cout << count;
    
    return NULL;

} 

*/
