#ifndef BINTREE_H
#define BINTREE_H
#include<stdio.h>
#include<vector>
#include<iostream>

using namespace std;

class BinTree
{
	public:
	//BinTree():tr(NULL),level(NULL), numNodes(0){}
	BinTree(int n):numNodes(0){Create(n);}
	~BinTree();
	std::vector<int> GetChildren(int ID);
	int GetSibling(int ID);
	/*bool IsLeaf(int ID);*/
	int GetLeftMostChild(int ID);
	int GetRightMostChild(int ID);
	std::vector<int> GetDescendents(int ID);
	int GetNumNodes() { return numNodes;}
	std::vector<int> GetLeaves() { return leaves;}
#ifdef DEBUG
	void Print();
#endif

	/* is a list to store the parent IDs. tr[0] contains the parent ID of node number 1, 
	 * tr[1] contains the parent ID of node number 2, and so on. The root node is represented by
	 * the ID, n, where n+1 is a perfect power of two (because it is a complete binary tree).
	 * tr[n-1] = 0 i.e. the parent ID of root node is 0 */
	int* tr; 

	/* is a list to store the level of each node starting from level 0 (root node) */
	int* level;

	/* is a list to store the IDs of leaf nodes */
	std::vector<int> leaves;

	int maxLevel; //the maximum level of a node in the tree for a given value of numNodes 
	int numNodes; //number of nodes in the complete binary tree.

	//helper functions to create the tree
	int BinaryPostorderTree(int parentID, int step, bool leftChild);
	void Create(int n);
};

#endif
