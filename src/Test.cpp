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

//define the value that triggers fmm acceleration. TODO: add detailed comment about for what values the acceleration is expected.
int fmmTrigger=17000; 
const char* testFile="sparseOut.txt"; // Pritesh: Band5.txt is there in Pnk system inside /hsseigen/ directory

template<typename T>
void PrintArray(T *Arr, int row, int col, const char* filename="output.txt")
{
	std::ofstream txtOut;
    txtOut.open(filename, std::ofstream::out | std::ofstream::app);
    for(unsigned int i = 0; i <(unsigned)row; i++)
    {
        for(int j=0; j < col; j++)
        {
            txtOut<<std::setprecision(12)<<Arr[j+i*col]<<"\n";
        }
    }
    txtOut.close();
}



int main(int argc, char* argv[])
{
	int n=16;
	int r=4;
	int MorB = 2; // to select which routine to use Band2HSS (2) or Mat2HSSsymm (1).
	int w = 4;	//  band of the matrix, can be changes according to the requirement.
	int nProc = 1; // no of procs to use
	if(argc!=7){
		printf("Usage: ./Test <filename> <matrix_size> <diagblock_size> <Band2HSS or Mat2HSSsym> <bandwidth> <no of processor>\n");
		
	}else{
		testFile=argv[1];
		n=atoi(argv[2]);
		r=atoi(argv[3]);
		MorB = atoi(argv[4]);
		w = atoi(argv[5]);
		nProc = atoi(argv[6]);
	}

	BinTree* bt=NULL;
	double * A=NULL;
	//create a banded matrix
	int status = MakeBand(n,w,&A);
	if(status) 
		exit(-1);
	
	
	int numNodes = 0;
	

	/*(you can either reuse the tree created earlier or let the call to NPart create a new tree based on the size of the partition specified.
	 * Arguments of NPart: n is the number of rows/columns in an input matrix. r is the number of rows in a partition (horizontal) of the matrix 
	 * Number of leaves = n/r. Num nodes in the tree = num leaves* 2 - 1*/
	int *m=NULL;
	int mSize;
	NPart(n, r, &bt, &m, mSize, numNodes);



	GEN *hss = HssGenerators(A, n*n, bt, m , mSize, w, MorB);
	
	
	SDC* res = superDC(hss, bt, m, mSize, nProc);

	//computing matvec using res->Q
	//x = zeros(n,1);
	//x(2) = 1;
	//y = superdcmv(Q,x,0, N);
	int k = bt->GetNumNodes(); //k is the ID of the root node.
	double* x=new double[32];//double[(res->qSizes[k-1]).second];
	//for(int i=0;i<(res->qSizes[k-1]).second;i++)
	for(int i=0;i<32;i++)
		x[i]=0;
	x[0]=1;

        int columns = 32;//(res->qSizes[k-1]).second;
        std::pair<int,int> xSize = {columns,1};
	//the index parameter (k) of superdcmv is expected to be the ID of the node. Since the functions called within superdcmv and their callees use the index to offset into qSize and Q, the index is decremented by 1.
	double* matvecprod=superdcmv(res->Q, res->qSizes, x, &xSize, m, mSize, bt, k-1, 0, fmmTrigger);
	PrintArray<double>(matvecprod, columns, 1, "eigvec.txt"); 

	delete[] m;
	
	return 0;
}
