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
        cout << "Using Band2HSS\n";
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

    int numD=N/2;
    int numURB=N-1;
    
    init->dSizes = new std::pair<int,int>[numD];
    init->uSizes = new std::pair<int,int>[numURB];
    init->rSizes = new std::pair<int,int>[numURB];
    init->bSizes = new std::pair<int,int>[numURB];

    init->D = new double*[numD];
    init->U = new double*[numURB];
    init->R = new double*[numURB];
    init->B = new double*[numURB];

    return init;
}

#ifdef DIST
#include <mpi.h>
void sendGen(GEN* send, MPI_Comm process_grid, int N){
    int numD = N/2;
    int numURB = N-1;

    int myrank=0;
    MPI_Comm_rank(process_grid, &myrank);
    // send Dsizes
    double *dszbuff = new double[numD*2];
    if (myrank==0){
        for (int i = 0; i < numD; i++)
        {
            dszbuff[i*2] = send->dSizes[i].first;
            dszbuff[i*2+1] = send->dSizes[i].second;
        }
    }
    MPI_Bcast(dszbuff, numD*2, MPI_DOUBLE, 0, process_grid);
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
        if (myrank!=0)
        {
            send->D[i] = new double[send->dSizes[i].first*send->dSizes[i].second];
        }
        MPI_Bcast(send->D[i], send->dSizes[i].first*send->dSizes[i].second, MPI_DOUBLE, 0, process_grid);
    }


    // Send U sizes
    double *Uszbuff = new double[numURB*2];
    if (myrank==0){
        for (int i = 0; i < numURB; i++)
        {
            Uszbuff[i*2] = send->uSizes[i].first;
            Uszbuff[i*2+1] = send->uSizes[i].second;
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
        if (myrank!=0)
        {
            send->U[i] = new double[send->uSizes[i].first*send->uSizes[i].second];
        }
        MPI_Bcast(send->U[i], send->uSizes[i].first*send->uSizes[i].second, MPI_DOUBLE, 0, process_grid);
    }
    


    // Send R sizes
    double *Rszbuff = new double[numURB*2];
    if (myrank==0){
        for (int i = 0; i < numURB; i++)
        {
            Rszbuff[i*2] = send->rSizes[i].first;
            Rszbuff[i*2+1] = send->rSizes[i].second;
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
        if (myrank!=0)
        {
            send->R[i] = new double[send->rSizes[i].first*send->rSizes[i].second];
        }
        MPI_Bcast(send->R[i], send->rSizes[i].first*send->rSizes[i].second, MPI_DOUBLE, 0, process_grid);
    }

    // Send B sizes
    double *Bszbuff = new double[numURB*2];
    if (myrank==0){
        for (int i = 0; i < numURB; i++)
        {
            Bszbuff[i*2] = send->bSizes[i].first;
            Bszbuff[i*2+1] = send->bSizes[i].second;
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
        if (myrank!=0)
        {
            send->B[i] = new double[send->bSizes[i].first*send->bSizes[i].second];
        }
        MPI_Bcast(send->B[i], send->bSizes[i].first*send->bSizes[i].second, MPI_DOUBLE, 0, process_grid);
    }
    
    
}
#endif