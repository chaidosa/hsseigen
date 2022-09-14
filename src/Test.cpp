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

//define the value that triggers fmm acceleration. TODO: add detailed comment about for what values the acceleration is expected.
int fmmTrigger=17000; 
const char* testFile="sparseOut.txt"; // Pritesh: Band5.txt is there in Pnk system inside /hsseigen/ directory

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
	int k = bt->GetNumNodes();
	double* x=new double[(res->qSizes[k-1]).second];
	for(int i=0;i<(res->qSizes[k-1]).second;i++)
		x[i]=0;
	x[0]=1;
	superdcmv(res->Q, res->qSizes, x, bt, 1, fmmTrigger);
	delete[] m;
	
	return 0;
}
