#include<iostream>
#include<stdio.h>
#include"makeband.h"
#include"BinTree.h"
#include"NPart.h"
#include"test_mat2hsssym.h"
#include "QR.h"
#include "superDC.h"
#include "band2hss.h"
#include <sys/time.h>
#include<assert.h>

char* testFile="sparseOut.txt";

int main(int argc, char* argv[])
{
	int n=126, band_width=2;
	int r=32;
	int MorB = 1; // to select which routine to use Band2HSS (2) or Mat2HSSsymm (1).
	int w = 0;	//  band of the matrix, can be changes according to the requirement.
	int nProc = 1;
	if(argc!=7){
		printf("Usage: ./Test <filename> <matrix_size> <diagblock_size> <Band2HSS or Mat2HSSsym> <bandwidth>\n");
		printf("Default values will be used: filename=sparseOut.txt, n=126, diagblock_size=32, Mat2hsssym, 0\n");
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
	int status = MakeBand(n,band_width,&A);
	if(status) 
		exit(-1);
	
	
	int numNodes = 0;
	
#ifdef DEBUG
	//Testing makeband
	cout<<"========Input Matrix================"<<endl;
	if(A)
	{
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				printf("%1.2f ",A[j+i*n]);
			}
			printf("\n");
		}
	}
	cout<<"========Input Matrix================"<<endl;
	//Testing BinTree
	numNodes = 15;
	bt = new BinTree(numNodes);
	if(bt->GetNumNodes() == 0)
		exit(-1);
	bt->Print();
	cout<<"Sibling of 7: "<<bt->GetSibling(7)<<endl;
	cout<<"Sibling of 14: "<<bt->GetSibling(14)<<endl;
	cout<<"Sibling of 6: "<<bt->GetSibling(6)<<endl;
	cout<<"Sibling of 4: "<<bt->GetSibling(4)<<endl;
	cout<<"Sibling of 1: "<<bt->GetSibling(1)<<endl;
	
	cout<<"========================"<<endl;
	std::vector<int> children = bt->GetChildren(7);
	if(children.size() == 2)
		cout<<"children of 7: "<<children[0]<<", "<<children[1]<<endl;
	else
		cout<<"No children found for 7"<<endl;

	children.clear();
	children = bt->GetChildren(14);
	if(children.size() == 2)
		cout<<"children of 14: "<<children[0]<<", "<<children[1]<<endl;
	else
		cout<<"No children found for 14"<<endl;

	children.clear();
	children = bt->GetChildren(6);
	if(children.size() == 2)
		cout<<"children of 6: "<<children[0]<<", "<<children[1]<<endl;
	else
		cout<<"No children found for 6"<<endl;

	children.clear();
	children = bt->GetChildren(5);
	if(children.size() == 2)
		cout<<"children of 5: "<<children[0]<<", "<<children[1]<<endl;
	else
		cout<<"No children found for 5"<<endl;

	children.clear();
	children = bt->GetChildren(1);
	if(children.size() == 2)
		cout<<"children of 1: "<<children[0]<<", "<<children[1]<<endl;
	else
		cout<<"No children found for 1"<<endl;
	
	cout<<"========================"<<endl;
	cout<<"Left most child of 15: "<<bt->GetLeftMostChild(15)<<endl;
	cout<<"Left most child of 14: "<<bt->GetLeftMostChild(14)<<endl;
	cout<<"Left most child of 7: "<<bt->GetLeftMostChild(7)<<endl;
	cout<<"Left most child of 10: "<<bt->GetLeftMostChild(10)<<endl;
	cout<<"Left most child of 2: "<<bt->GetLeftMostChild(2)<<endl;
	cout<<"========================"<<endl;
	
	cout<<"Right most child of 15: "<<bt->GetRightMostChild(15)<<endl;
	cout<<"Right most child of 14: "<<bt->GetRightMostChild(14)<<endl;
	cout<<"Right most child of 7: "<<bt->GetRightMostChild(7)<<endl;
	cout<<"Right most child of 10: "<<bt->GetRightMostChild(10)<<endl;
	cout<<"Right most child of 2: "<<bt->GetRightMostChild(2)<<endl;
	cout<<"========================"<<endl;

	std::vector<int> desc = bt->GetDescendents(15);
	if(desc.size() > 0) {
		cout<<"descendants of 15: ";
		for(int i=0;i<desc.size();i++) {
			cout<<desc[i];
			if(i != desc.size()-1)
				cout<<" ";
			else
				cout<<endl;
		}
	}
	else
		cout<<"No descendants found for 15"<<endl;
	
	desc.clear();
	desc = bt->GetDescendents(14);
	if(desc.size() > 0) {
		cout<<"descendants of 14: ";
		for(int i=0;i<desc.size();i++) {
			cout<<desc[i];
			if(i != desc.size()-1)
				cout<<" ";
			else
				cout<<endl;
		}
	}
	else
		cout<<"No descendants found for 14"<<endl;

	desc.clear();
	desc = bt->GetDescendents(6);
	if(desc.size() > 0) {
		cout<<"descendants of 6: ";
		for(int i=0;i<desc.size();i++) {
			cout<<desc[i];
			if(i != desc.size()-1)
				cout<<" ";
			else
				cout<<endl;
		}
	}
	else
		cout<<"No descendants found for 6"<<endl;

	desc.clear();
	desc = bt->GetDescendents(7);
	if(desc.size() > 0) {
		cout<<"descendants of 7: ";
            //std::copy(E, E+resDvd->dSizes[i].first, ee);
		for(int i=0;i<desc.size();i++) {
			cout<<desc[i];
			if(i != desc.size()-1)
				cout<<" ";
			else
				cout<<endl;
		}
	}
	else
		cout<<"No descendants found for 7"<<endl;
	
	desc.clear();
	desc = bt->GetDescendents(2);
	if(desc.size() > 0) {
		cout<<"descendants of 2: ";
		for(int i=0;i<desc.size();i++) {
			cout<<desc[i];
			if(i != desc.size()-1)
				cout<<" ";
			else
				cout<<endl;
		}
	}
	else
		cout<<"No descendants found for 2"<<endl;
	cout<<"========================"<<endl;
#endif
	/*(you can either reuse the tree created earlier or let the call to NPart create a new tree based on the size of the partition specified.
	 * Arguments of NPart: n is the number of rows/columns in an input matrix. r is the number of rows in a partition (horizontal) of the matrix 
	 * Number of leaves = n/r. Num nodes in the tree = num leaves* 2 - 1*/
	int *m=NULL;
	int mSize;
	NPart(n, r, &bt, &m, mSize, numNodes);

#if DEBUG	
	//Testing NPart
	if(bt && m)
	{
		if(numNodes == 0)
			bt->Print();
		cout<<"========partition sizes==============="<<endl;
		for(int i=0;i<mSize;i++) {
			cout<<m[i];
			if(i!=mSize-1)
				cout<<" ";
			else
				cout<<endl;
		}
		cout<<"========partition sizes==============="<<endl;
	}
#endif

	GEN *hss = HssGenerators(A, n*n, bt, m , mSize, w, MorB);
	
	
	SDC* res = superDC(hss, bt, m, mSize, nProc);
	//res->~SDC();
	delete[] m;
	
	return 0;
}
