#ifndef BINTREE_H
#define BINTREE_H
#include<stdio.h>
#include<vector>
#include<iostream>

using namespace std;

class BinTree
{
	public:
	
	BinTree(int n):numNodes(0){Create(n);}
	//~BinTree();
	std::vector<int> GetChildren(int ID);

	std::vector<int> GetTreeDesc();
	int GetNumNodes();

	vector<int> tr; 

	int numNodes; //number of nodes in the complete binary tree.

	//helper functions to create the tree
	void Create(int n);
};

#endif
