#include<string.h>
#include<math.h>
#include<algorithm>
#include<assert.h>
#include "BinTree.h"
#include "band2hss.h"
#include "Generators.h"
#include "iostream"
#include "test_mat2hsssym.h"

using namespace std;


GEN* HssGenerators(double *A, int aSize, BinTree* bt, int* m, int mSize, int w, int MorB){
    
    GEN *res;
    if(MorB == 2){
        // cout << "Using Band2HSS\n";
        res = band2hss(A, aSize, bt, m, mSize, w);
    }
    else if(MorB == 1){
        cout << "Using Mat2Hssym\n";
        res = t_mat2hsssym(A, aSize, bt, m, mSize);      
    }
    else{
        cout<<"Something went wrong!!!, Try again\n";
        assert(false);
    }


    return res;
}

GEN* InitGen(int N){
    GEN* init = new GEN();

    int numD=N;
    int numURB=N;
    int n=N;
    
    init->dSizes = new std::pair<int,int>[numD];
    for(int i = 0; i < n; i++)
		init->dSizes[i]=std::make_pair(0,0); 
    init->uSizes = new std::pair<int,int>[numURB];
    for(int i = 0; i < n; i++)
		init->uSizes[i]=std::make_pair(0,0); 
    init->rSizes = new std::pair<int,int>[numURB];
    for(int i = 0; i < n; i++)
		init->rSizes[i]=std::make_pair(0,0); 
    init->bSizes = new std::pair<int,int>[numURB];
    for(int i = 0; i < n; i++)
		init->bSizes[i]=std::make_pair(0,0); 

    init->D = new double*[numD];
    for (int i = 0; i < n; i++) init->D[i] = NULL;
    init->U = new double*[numURB];
    for (int i = 0; i < n; i++) init->U[i] = NULL;
    init->R = new double*[numURB];
    for (int i = 0; i < n; i++) init->R[i] = NULL;
    init->B = new double*[numURB];
    for (int i = 0; i < n; i++) init->B[i] = NULL;

    return init;
}

#if defined(DIST) || defined(HYBRD)
#include <mpi.h>
void sendGen(GEN* send, MPI_Comm process_grid, int N){
    int numD = N;
    int numURB = N;


    int myrank=0;
    MPI_Comm_rank(process_grid, &myrank);
    // cout << "Rank is " << myrank << "\n";
    // send Dsizes
    int *dszbuff = new int[numD*2];
    if (myrank==0){
        for (int i = 0; i < numD; i++)
        {
            if (send->D[i]==NULL)
            {
                dszbuff[i*2] = -1;
                dszbuff[i*2+1] = -1;
            } else {
                dszbuff[i*2] = send->dSizes[i].first;
                dszbuff[i*2+1] = send->dSizes[i].second;
            }
        }
    }
    MPI_Bcast(dszbuff, numD*2, MPI_INT, 0, process_grid);
    if (myrank!=0)
    {
        for (int i = 0; i < numD; i++)
        {
            send->dSizes[i].first = dszbuff[i*2];
            send->dSizes[i].second = dszbuff[i*2+1];  
        }
    }
    delete[] dszbuff;
    
    // send D
    for (int i = 0; i < numD; i++)
    {   
        int recv_elements = send->dSizes[i].first*send->dSizes[i].second;
        if (myrank!=0)
        {   
            if (send->dSizes[i].first == -1 || send->dSizes[i].second == -1)
            {
                send->D[i] = NULL;
                send->dSizes[i].first=0;
                send->dSizes[i].second=0;
                recv_elements=0;
            } else {    
                send->D[i] = new double[send->dSizes[i].first*send->dSizes[i].second];
                recv_elements = send->dSizes[i].first*send->dSizes[i].second;
            }
        }

        MPI_Bcast(send->D[i], recv_elements, MPI_DOUBLE, 0, process_grid);
    }

    // Send U sizes
    double *Uszbuff = new double[numURB*2];
    if (myrank==0){
        for (int i = 0; i < numURB; i++)
        {   
            if (send->U[i]==NULL)
            {
                Uszbuff[i*2] = -1;
                Uszbuff[i*2+1] = -1;
            } else {
                Uszbuff[i*2] = send->uSizes[i].first;
                Uszbuff[i*2+1] = send->uSizes[i].second;
            }
        }
    }
    MPI_Bcast(Uszbuff, numURB*2, MPI_DOUBLE, 0, process_grid);
    if (myrank!=0)
    {
        for (int i = 0; i < numURB; i++)
        {
            send->uSizes[i].first = Uszbuff[i*2];
            send->uSizes[i].second = Uszbuff[i*2+1];  
        }
    }
    delete[] Uszbuff;

    // Send U
    for (int i = 0; i < numURB; i++)
    {   
        int recv_elements = send->uSizes[i].first*send->uSizes[i].second;
        if (myrank!=0)
        {   
            if (send->uSizes[i].first == -1 || send->uSizes[i].second == -1)
            {
                send->U[i] = NULL;
                send->uSizes[i].first=0;
                send->uSizes[i].second=0;
                recv_elements = 0;
            } else{
                send->U[i] = new double[send->uSizes[i].first*send->uSizes[i].second];
                recv_elements = send->uSizes[i].first*send->uSizes[i].second;
            }
        }

        MPI_Bcast(send->U[i], recv_elements, MPI_DOUBLE, 0, process_grid);
    }
    


    // Send R sizes
    double *Rszbuff = new double[numURB*2];
    if (myrank==0){
        for (int i = 0; i < numURB; i++)
        {   
            if (send->R[i]==NULL)
            {
                Rszbuff[i*2] = -1;
                Rszbuff[i*2+1] = -1;
            } else {
                Rszbuff[i*2] = send->rSizes[i].first;
                Rszbuff[i*2+1] = send->rSizes[i].second;
            }
        }
    }
    MPI_Bcast(Rszbuff, numURB*2, MPI_DOUBLE, 0, process_grid);
    if (myrank!=0)
    {
        for (int i = 0; i < numURB; i++)
        {
            send->rSizes[i].first = Rszbuff[i*2];
            send->rSizes[i].second = Rszbuff[i*2+1];  
        }
    }
    delete[] Rszbuff;

    // Send R
    for (int i = 0; i < numURB; i++)
    {   
        int recv_elements = send->rSizes[i].first*send->rSizes[i].second;
        if (myrank!=0)
        {   
            if (send->rSizes[i].first == -1 || send->rSizes[i].second == -1)
            {
                send->R[i] = NULL;
                send->rSizes[i].first=0;
                send->rSizes[i].second=0;
                recv_elements=0;
            } else {
                send->R[i] = new double[send->rSizes[i].first*send->rSizes[i].second];
                recv_elements = send->rSizes[i].first*send->rSizes[i].second;
            }
        }

        MPI_Bcast(send->R[i], recv_elements, MPI_DOUBLE, 0, process_grid);
    }

    // Send B sizes
    double *Bszbuff = new double[numURB*2];
    if (myrank==0){
        for (int i = 0; i < numURB; i++)
        {   
            if (send->B[i]==NULL)
            {
                Bszbuff[i*2] = -1;
                Bszbuff[i*2+1] = -1;
            } else {
                Bszbuff[i*2] = send->bSizes[i].first;
                Bszbuff[i*2+1] = send->bSizes[i].second;
            }
        }
    }
    MPI_Bcast(Bszbuff, numURB*2, MPI_DOUBLE, 0, process_grid);
    if (myrank!=0)
    {
        for (int i = 0; i < numURB; i++)
        {
            send->bSizes[i].first = Bszbuff[i*2];
            send->bSizes[i].second = Bszbuff[i*2+1];  
        }
    }
    delete[] Bszbuff;

    // send B
    for (int i = 0; i < numURB; i++)
    {
        int recv_elements = send->bSizes[i].first*send->bSizes[i].second;
        if (myrank!=0)
        {   
            if (send->bSizes[i].first == -1 || send->bSizes[i].second == -1)
            {
                send->B[i] = NULL;
                send->bSizes[i].first=0;
                send->bSizes[i].second=0;
                recv_elements=0;
            } else {
                send->B[i] = new double[send->bSizes[i].first*send->bSizes[i].second];
                recv_elements = send->bSizes[i].first*send->bSizes[i].second;
            }
        }

        MPI_Bcast(send->B[i], recv_elements, MPI_DOUBLE, 0, process_grid);
    }
    
    // cout << "Completed B\n";    
}
#endif