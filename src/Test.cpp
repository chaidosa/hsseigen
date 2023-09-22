#include<iostream>
#include<stdio.h>
#include"makeband.h"
#include"BinTree.h"
#include"NPart.h"
#include"test_mat2hsssym.h"
#include "QR.h"
#include "superDC.h"
#include "band2hss.h"
#include"superdcmv.h"
#include <sys/time.h>
#include<assert.h>
#include<fstream>
#include<iomanip>
#include<numeric>
#include<algorithm>

//define the value that triggers fmm acceleration. TODO: add detailed comment about for what values the acceleration is expected.
int fmmTrigger=17000; 
const char* testFile="sparseOut.txt"; // Pritesh: Band5.txt is there in Pnk system inside /hsseigen/ directory

template<typename T>
void PrintArray(T *Arr, int row, int col, const char* filename="output.txt")
{
	std::ofstream txtOut(filename);
    for(unsigned int i = 0; i <(unsigned)row; i++)
    {
        for(int j=0; j < col; j++)
        {
            txtOut<<std::setprecision(12)<<Arr[j+i*col]<<"\n";
        }
    }
    txtOut.close();
}


#if defined(DIST) || defined(HYBRD)
// extern "C"{
#include <mpi.h>
#include <cmath>
#include <cstring>
// }

void fill_arr(double *A, int n){

    //srand( (unsigned)time(NULL) );

    for (int i = 0; i < n; i++)

    {

        A[i] = (double) rand()/10e+8;

    }
}
#endif

int main(int argc, char* argv[])
{

#if defined(DIST) || defined(HYBRD)
	MPI_Init( &argc, &argv);
#endif



	int n=16;
	int r=4;
	int MorB = 2; // to select which routine to use Band2HSS (2) or Mat2HSSsymm (1).
	int w = 4;	//  band of the matrix, can be changes according to the requirement.
	int nProc = 1; // no of procs to use
	if(argc!=7){
		printf("Usage: ./Test <filename> <matrix_size> <diagblock_size> <Band2HSS or Mat2HSSsym> <bandwidth> <no of processor>\n");
		// testFile = "/home/c200010021/hssdata/halfb5_32k.txt";
		// n = 32768;
		// r = 64;
		// MorB = 2;
		// w = 5;
		// nProc = 1;
	}else{
		testFile=argv[1]; 		// filename
		n=atoi(argv[2]); 		// matrix size
		r=atoi(argv[3]); 		// diag_block size
		MorB = atoi(argv[4]);   // band2hss 
		w = atoi(argv[5]);		// bandwidth
		nProc = atoi(argv[6]);	// num. Processors
	}

	int myrank = 0;
	#if defined(DIST) || defined(HYBRD)
		int procs;
    	MPI_Comm_size(MPI_COMM_WORLD, &procs);
		// r = (int)n/procs;
	    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	#endif

	#ifdef HYBRD
		r = atoi(argv[3]);
	#endif


	BinTree* bt=NULL;
	double * A=NULL; // The matrix A from test_input.txt of size n*n

// create a banded matrix
#if defined(DIST)
	int numElemsInaRow = n;
	int size = (w + w + 1) * numElemsInaRow;
	A = new double[size];
	memset(A, 0, sizeof(double) * size);

	for (int i = 0; i < (w + 1); i++)
	{
		if (i != w)
		{
			int nn = numElemsInaRow - w + i;
			double *temp_arr = new double[nn];
			fill_arr(temp_arr, nn);

			double *Point_a = A + (i * numElemsInaRow) + (w - i);
			double *Point_b = A + (2 * w - i) * numElemsInaRow;

			memcpy(Point_a, temp_arr, sizeof(double) * nn);
			memcpy(Point_b, temp_arr, sizeof(double) * nn);
			delete[] temp_arr;
		}

		else
		{

			double *temp_arr = new double[numElemsInaRow];
			fill_arr(temp_arr, numElemsInaRow);
			double *point = A + i * numElemsInaRow;

			memcpy(point, temp_arr, sizeof(double) * numElemsInaRow);
			delete[] temp_arr;
		}
	}
#else
	int status = MakeBand(n, w, &A); // Makes band matrix of n rows and w bandwidth
	if (status)
		exit(-1);
#endif

	int numNodes = 0;

	/*(you can either reuse the tree created earlier or let the call to NPart create a new tree based on the size of the partition specified.
	 * Arguments of NPart: n is the number of rows/columns in an input matrix. r is the number of rows in a partition (horizontal) of the matrix 
	 * Number of leaves = n/r. Num nodes in the tree = num leaves* 2 - 1*/
	int *m=NULL;
	int mSize;
	NPart(n, r, &bt, &m, mSize, numNodes);

	GEN *hss;
	#if defined(DIST) || defined(HYBRD)
	MPI_Barrier(MPI_COMM_WORLD);

	if (myrank==0){
		hss = HssGenerators(A, n*n, bt, m , mSize, w, MorB); 
	} else {
		hss = InitGen(bt->GetNumNodes()); 
	}
	
	// cout << "\nStarted Hssgen\n";
	sendGen(hss, MPI_COMM_WORLD, bt->GetNumNodes()); // check
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "\nCompletes Hssgen\n";
	#else
	cout << "Completed Hss tree\n";
	hss =  HssGenerators(A, n*n, bt, m , mSize, w, MorB);
	#endif

	SDC* res; 

	#if defined(DIST) || defined(HYBRD)
	/**
     * Create a 2d grid of processors PxQ where LCM(P,Q) = 1
     */
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    // int p = 2, q = 3;

	// if (p*q != nprocs)
	// {
	// 	cout << "Error PxQ must be equal to nprocs" << "\n";
	// 	MPI_Abort(MPI_COMM_WORLD, 0);
	// }
	

	// MPI_Comm process_grid;
    // int dims[2] = {p, q}; // integer array of size ndims specifying the number of processes in each dimension
    // int periods[2] = {0, 0};

    // MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &process_grid);

	// cout << "Reached dsuperdc" << "\n";
	#endif

	#ifdef DIST
	res = dsuperDc(hss, bt, m, mSize, nProc, MPI_COMM_WORLD);
	#endif

	#ifdef HYBRD
	res = HybridSuperDC(hss, bt, m, mSize, nProc, MPI_COMM_WORLD);
	#endif

	#if defined(DIST) || defined(HYBRD)
	if (myrank==0){
	#else
	res = superDC(hss, bt, m, mSize, nProc);
	#endif
	//this sorting of eigenvalues from smallest to largest is necessary and the index values
	vector<int> I(res->lSize);
 	std::iota(I.begin(),I.end(),0); 
	std::sort( I.begin(),I.end(), [&](int i,int j){return res->L[i]<res->L[j];} );
	//computing matvec using res->Q
	//x = zeros(n,1);
	//x(2) = 1;
	//y = superdcmv(Q,x,0, N);
	int k = bt->GetNumNodes(); //k is the ID of the root node.
	double* x=new double[n];//double[(res->qSizes[k-1]).second];
	//for(int i=0;i<(res->qSizes[k-1]).second;i++)
	for(int i=0;i<n;i++)
		x[i]=0;
	x[0]=1;

        int columns = n;//(res->qSizes[k-1]).second;
        std::pair<int,int> xSize = {columns,1};
	assert(res->lSize==n);
	//the index parameter (k) of superdcmv is expected to be the ID of the node. Since the functions called within superdcmv and their callees use the index to offset into qSize and Q, the index is decremented by 1.
	double* matvecprod=superdcmv(res->Q, res->qSizes, x, &xSize, m, mSize, I.data(), bt, k-1, 0, fmmTrigger);
	PrintArray<double>(matvecprod, columns, 1, "eigvec.txt"); 

	delete[] m;

	#if	defined(DIST) || defined(HYBRD)
		}
		MPI_Finalize();
	#endif

	return 0;
}
