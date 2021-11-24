#include<math.h>
#include<assert.h>
#include"BinTree.h"
void BinTree::Create(int n)
{
	// Not a complete binary tree;
	if((n+1) & n) {
		cout<<"Error cannot create a complete binary tree with "<<n<<"nodes"<<endl;
		return;
	}
	numNodes=n;
	cout<<"creating tree with "<<numNodes<<"nodes"<<endl;

	//reserve space to hold parent IDs of n vertices
	tr = new int[numNodes];

	//reserve space to hold level numbers of n vertices
	level = new int[numNodes];

	//set parent ID of root node
	tr[numNodes-1]=0;

	//set level of root node
	level[numNodes-1]=0;

	//step is one more than the number of nodes in a subtree rooted at left/right child
	int step = (numNodes+1)/2;

	//maximum possible level (levels starting from 0. Hence, -1.
	maxLevel = log2(numNodes+1) - 1;

	//create IDs of left subtree
	int leftID = BinaryPostorderTree(numNodes, step, true);

	//create IDs of right subtree
	int rightID = BinaryPostorderTree(numNodes, step, false);

	//make root node the parent of left and right children
	tr[rightID-1]=numNodes;
	tr[leftID-1]=numNodes;

	return;
}

BinTree::~BinTree()
{
	if(tr)
		delete [] tr;
	if(level)
		delete [] level;
}

int BinTree::BinaryPostorderTree(int parentID, int step, bool leftChild)
{
	int myID;
	if(!leftChild)
		myID = parentID-1;
	else
		myID = parentID - step;

	level[myID-1] = level[parentID-1] + 1;
	if(step == 2)
	{
		leaves.push_back(myID);
		return myID;
	}

	step=step/2;
	int leftID = BinaryPostorderTree(myID, step, true);
	int rightID = BinaryPostorderTree(myID, step, false);
	tr[rightID-1] = myID;
	tr[leftID-1] = myID;
	return myID;
}

std::vector<int> BinTree::GetChildren(int ID)
{
	std::vector<int> ret;
	for(int i=0;i<numNodes;i++)
	{
		if(tr[i] == ID)
		{
			ret.push_back(i+1);
			if(ret.size() == 2)
				break;
		}
	}
	return ret;
}

int BinTree::GetSibling(int ID)
{
	int sib=-1;
	if((ID >= numNodes) || (ID <= 0))
		return sib;

	int parentID = tr[ID-1];
	int step = (numNodes+1)/(1<<(level[ID-1]));
	if(ID == parentID-1)
	{
		//I am rightChild. Get leftSibling ID.
		sib = parentID-step;
	}
	else
	{
		//I am leftChild. Get rightSibling ID.
		sib = parentID-1;
	}

	return sib;
}

/*bool BinTree::IsLeaf(int ID)
{
	std::vector<int>::iterator iter = std::find(leaves.begin(), leaves.end(), ID);
	if(iter == leaves.end())
		return false;
	else
		return true;
}*/

int BinTree::GetLeftMostChild(int ID)
{
	int curLevel = level[ID-1];
	int numNodesInSubtreeOfChild = 1<<(maxLevel-curLevel);

	while(level[ID-1] < maxLevel) {
		ID = ID-numNodesInSubtreeOfChild; //since leftChildID is equal to parentID-step.
		numNodesInSubtreeOfChild /= 2;
	}
	return ID;
}

int BinTree::GetRightMostChild(int ID)
{
	while(level[ID-1] < maxLevel)
		ID = ID-1; //since the rightChildID is equal to parentID-1.
	return ID;
}

std::vector<int> BinTree::GetDescendents(int ID)
{
	int curLevel = level[ID-1];
	//Number of descendents = totalNodesInSubtree-1.
	int totalNodesInSubtree = (1<<(maxLevel-curLevel+1))-1;
	std::vector<int> ret(totalNodesInSubtree-1);
	int leftMostID = GetLeftMostChild(ID);
	int numAdded=0;
	assert(ID-leftMostID == ret.size());
	for(int i=leftMostID;i<ID;i++)
		ret[numAdded++]=i;
	return ret;
}

std::vector<int> BinTree::GetTreeDesc()
{
	std::vector<int>ret(numNodes+1);
	ret[0] = 0;
	for(int i=1;i<=numNodes;i++){
		ret[i] = GetLeftMostChild(i);
	}
	return ret;
}

#ifdef DEBUG
void BinTree::Print() {
	cout<<"numNodes: "<<numNodes<<" maxLevel: "<<maxLevel<<endl;
	cout<<"----------level--------"<<endl;
	for(int i=0;i<numNodes;i++) {
		cout<<level[i];
		if( i!= numNodes - 1)
			cout<<" ";
		else
			cout<<endl;
	}
	cout<<"----------level--------"<<endl;
	cout<<"----------parent IDs--------"<<endl;
	for(int i=0;i<numNodes;i++) {
		cout<<tr[i];
		if( i!= numNodes - 1)
			cout<<" ";
		else
			cout<<endl;
	}
	cout<<"----------parent IDs--------"<<endl;
	cout<<"----------Leaf Nodes--------"<<endl;
	for(int i=0;i<leaves.size();i++) {
		cout<<leaves[i];
		if( i!= leaves.size() - 1)
			cout<<" ";
		else
			cout<<endl;
	}
	cout<<"----------Leaf Nodes--------"<<endl;
}
#endif
