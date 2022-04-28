#include<iostream>
#include<stdio.h>
#include"makeband.h"
#include"BinTree.h"
#include"NPart.h"
#include"test_mat2hsssym.h"
#include "QR.h"
#include "superDC.h"
#include "band2hss.h"
#include "Generators.h"
#include <sys/time.h>
#include<assert.h>
#include <fstream>

int main()
{
	int n=4096, band_width=2;
	int r=256;
	int MorB = 1; // to select which routine to use Band2HSS (2) or Mat2HSSsymm (1).
	int w = 0;	//  band of the matrix, can be changes according to the requirement.

	BinTree* bt=NULL;
	
	
	int numNodes = 0;
	
	/*(you can either reuse the tree created earlier or let the call to NPart create a new tree based on the size of the partition specified.
	 * Arguments of NPart: n is the number of rows/columns in an input matrix. r is the number of rows in a partition (horizontal) of the matrix 
	 * Number of leaves = n/r. Num nodes in the tree = num leaves* 2 - 1*/
	int *m=NULL;
	int mSize;
	NPart(n, r, &bt, &m, mSize, numNodes);
	int N = bt->numNodes;
	GEN *hss = new GEN();
	double **D = new double*[bt->numNodes];
	std::pair<int, int>*dSizes = new std::pair<int, int>[N];

	double **U = new double*[bt->numNodes];
	std::pair<int, int>*uSizes = new std::pair<int, int>[N];

	double **R = new double*[bt->numNodes];
	std::pair<int, int>*rSizes = new std::pair<int, int>[N];

	double **B = new double*[bt->numNodes];
	std::pair<int, int>*bSizes = new std::pair<int, int>[N];
	
	std::ifstream ifD ("genD4.txt", std::ifstream::in);
	if(!ifD.is_open()){
		cout << "Error in opening test file"<<endl;
		assert(false);
	}
	while(!ifD.eof()){
		int i, row, col;
		ifD >> i;
		ifD >> row;
		ifD >> col;
		dSizes[i-1] = make_pair(row, col);
		D[i-1] = new double[row * col];

		for(int r = 0; r < row; r++){
			for(int c = 0; c < col; c++){
				ifD >> D[i-1][c + r*col];
			}
		}	
	}
	
	 std::ifstream ifU ("genU4.txt", std::ifstream::in);
        if(!ifU.is_open()){
                cout << "Error in opening test file"<<endl;
                assert(false);
        }
        while(!ifU.eof()){
                int i, row, col;
                ifU >> i;
                ifU >> row;
                ifU >> col;
                uSizes[i-1] = make_pair(row, col);
                U[i-1] = new double[row * col];

                for(int r = 0; r < row; r++){
                        for(int c = 0; c < col; c++){
                                ifU >> U[i-1][c + r*col];
                        }
                }
        }
      	
	 std::ifstream ifR ("genR4.txt", std::ifstream::in);
        if(!ifR.is_open()){
                cout << "Error in opening test file"<<endl;
                assert(false);
        }
        while(!ifR.eof()){
                int i, row, col;
                ifR >> i;
                ifR >> row;
                ifR >> col;
                rSizes[i-1] = make_pair(row, col);
                R[i-1] = new double[row * col];

                for(int r = 0; r < row; r++){
                        for(int c = 0; c < col; c++){
                                ifR >> R[i-1][c + r*col];
                        }
                }
        }



	 std::ifstream ifB ("genB4.txt", std::ifstream::in);
        if(!ifB.is_open()){
                cout << "Error in opening test file"<<endl;
                assert(false);
        }
        while(!ifB.eof()){
                int i, row, col;
                ifB >> i;
                ifB >> row;
                ifB >> col;
                bSizes[i-1] = make_pair(row, col);
                B[i-1] = new double[row * col];

                for(int r = 0; r < row; r++){
                        for(int c = 0; c < col; c++){
                                ifB >> B[i-1][c + r*col];
                        }
                }
        }

	hss->D = D;
	hss->dSizes = dSizes;

	hss->B = B;
        hss->bSizes = bSizes;

	hss->R = R;
        hss->rSizes = rSizes;

	hss->U = U;
        hss->uSizes = uSizes;

	SDC* res = superDC(hss, bt, m, mSize, 1);
///	res->~SDC();
//	delete[] m;
	
	return 0;
}
