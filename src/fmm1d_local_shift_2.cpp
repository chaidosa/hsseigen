#include<cstdlib>
#include<cassert>
#include<math.h>
#include<set>
#include"fmm_types.h"
#define PARTITION_THRESHOLD 512 //number of points in leaf-node
#define TAU 0.6 //also called separation ratio or opening radius (in case of tree representing 2D/3D points) 
using namespace std;
typedef pair<double*, int*> TreeNodeData;

Vertex* rootNode = NULL;
vector<Vertex*> allLeaves;
double eta0;

const Vertex* GetRootNode() { return rootNode;}
void SetRootNode(Vertex* node) { rootNode=node;}

int compare(const void* x, const void* y) {
	const double arg1 = *static_cast<const double*>(x);
	const double arg2 = *static_cast<const double*>(y);
	if (arg1 < arg2) return -1;
	if (arg1 > arg2) return 1;
	return 0;
}

/* x and y contain elements at indices [xLb,xUb) [yLb,yUb). center is the absolute value of coordinate of the center of the interval. rad is the radius if the interval is bisected. */
Vertex* ConstructSpatialTree(double *x, double *y, int xLb, int xUb, int yLb, int yUb, double center, double rad, int el, int er, Vertex* parent) {

	double c = center - rad/2; //new center
	double d = rad/2; //new radius
	double mid;
	bool flag=false;
	int i, j, k, l, countXLeft=0, countXRight=0, countYLeft=0, countYRight=0;
	Vertex* leftChild=NULL, *rightChild=NULL, *node=NULL;

	if((xLb > xUb) || (yLb > yUb))
		return NULL;

	int sizex = xUb-xLb, sizey = yUb-yLb; //calculate number of elements in the interval
	if ((sizex <= MAX_POINTS_IN_CELL) || (sizey <= MAX_POINTS_IN_CELL)){
		//if this function is called with arguments such that the number of elements in an interval is below a threshold, then stop subdividing the interval.
		node = new Vertex();
		node->xLeft = xLb;
		node->xRight = xUb;
		node->yLeft = yLb;
		node->yRight = yUb;
		node->isLeaf = true;
		node->level = parent->level + 1;
		node->center = center; //note: not the new center.
		node->radius = rad; //note: not the new radius
		node->eta = eta0/rad;
		return node;
	}

	//el is the flag indicating if left boundary is to be included. get count of how many elements end up in the left interval (after division of the current interval) 
	if(el==1){
		for(i=xLb;i<xUb;i++) {
			if(((center-rad)<=x[i]) && (x[i]<=center))
				countXLeft++;
		}
		for(j=yLb;j<yUb;j++) {
			if(((center-rad)<=y[j]) && (y[j]<=center))
				countYLeft++;

		}
	} else if (el == 0) {
		for(i=xLb;i<xUb;i++) {
			if(((center-rad)<x[i]) && (x[i]<=center))
				countXLeft++;
		}
		for(j=yLb;j<yUb;j++){
			if(((center-rad)<y[j]) && (y[j]<=center))
				countYLeft++;
		}
	}

	//el is the flag indicating if left boundary is to be included. get count of how many elements end up in the right interval (after division of the current interval) 
	if(er==1) {
		for(k=xLb;k<xUb;k++)
			if((center<x[k]) && (x[k]<=(center+rad)))
				countXRight++;
		for(l=yLb;l<yUb;l++)
			if((center<y[l]) && (y[l]<=(center+rad)))
				countYRight++;
	} else if (er == 0) {
		for(k=xLb;k<xUb;k++)
			if(((center<x[k]) && (x[k]<(center+rad))))
				countXRight++;
		for(l=yLb;l<yUb;l++);
			if((center<y[l]) && (y[l]<(center+rad)))
				countYRight++;
	}
	
	if((countXLeft == 0) || (countYLeft == 0) || (countXRight == 0) || (countYRight == 0)){
		countXLeft=0; countXRight=0; countYLeft=0; countYRight=0;
		//x and y intervals are empty after the bisected point. Adjust center and rad.
		if(((xUb-xLb) > 0) && ((yUb-yLb) > 0))
			mid = (y[yUb]+x[xLb])/2;
		else if ((xUb-xLb) > 0)
			mid = (center+rad+x[xLb]) / 2;
		else if  ((yUb-yLb) > 0)
			mid = (center-rad+y[yUb]) / 2;
		c = (mid + (center - rad))/2;
		d = (mid - (center - rad))/2;
		if(el==1){
			for(i=xLb;i<xUb;i++)
				if(((center-rad)<=x[i]) && (x[i]<=mid))
					countXLeft++;
			for(j=yLb;j<yUb;j++)
				if(((center-rad)<=y[j]) && (y[j]<=mid))
					countYLeft++;
		} else if (el == 0) {
			for(i=xLb;i<xUb;i++)
				if(((center-rad)<x[i]) && (x[i]<=mid))
					countXLeft++;
			for(j=yLb;j<yUb;j++)
				if(((center-rad)<y[j]) && (y[j]<=mid))
					countYLeft++;
		}

		if(er==1) {
			for(k=xLb;k<=xUb;k++)
				if((mid<x[k]) && (x[k]<=(center+rad)))
					countXRight++;
			for(l=yLb;l<=yUb;l++)
				if((mid<y[l]) && (y[l]<=(center+rad)))
					countYRight++;
		} else if (er == 0) {
			for(k=xLb;k<=xUb;k++)
				if((mid<x[k]) && (x[k]<(center+rad)))
					countXRight++;
			for(l=yLb;l<=yUb;l++)
				if((mid<y[l]) && (y[l]<(center+rad)))
					countYRight++;
		}
		flag=true;
	}

	if((countXLeft > 0) || (countYLeft > 0)) {
		node = new Vertex();
		node->parent = parent;
		node->level = parent->level + 1;
		node->xLeft = xLb;
		node->xRight = xUb;
		node->yLeft = yLb;
		node->yRight = yUb;
		leftChild = ConstructSpatialTree(x,y,xLb,xLb+countXLeft,yLb,yLb+countYLeft, c, d, el, 1, node); 
	}
	else {
		assert(0); //should not come here.
	}
	
	c = center + rad/2; //if(!flag) then the value of c.
	if(flag) {
		c = (mid + (center + rad))/2; //overwrite the value of c if flag is true.
		d = ((center + rad) - mid)/2; 
	}
	
	assert(node);
	assert(xLb+countXLeft == xUb-countXRight);
	assert(yLb+countYLeft == yUb-countYRight);
	rightChild = ConstructSpatialTree(x,y,xLb+countXLeft,xUb,yLb+countYLeft,yUb, c, d, 0, er, node); 
	node->left = leftChild;
	node->right = rightChild;
	node->center = center;
	node->radius = rad;
	node->eta = eta0/rad;
	return node;
}

/*This is a 1D problem. So, when the input points are sorted and when the interval that they fall into is bisected,  the resulting two sub-intervals can be identified using the index of the left-most and right-most point falling in the respective sub-interval.*/
bool AreAdjacent(const Vertex* node1, const Vertex* node2) {
	bool ret = false;
	if((node1!=NULL) && (node2!=NULL)) {
		if((node1->xRight == node2->xLeft) || (node1->xLeft == node2->xRight) || (node1->yRight == node2->yLeft) || (node1->yLeft == node2->yRight))
			ret=true;
	}
	return ret;
}

void UpdateNeighbors(Vertex* node) {
	
	Vertex* parent = node->parent;
	for(int i=0;i<parent->neighbors.size();i++) {	
		Vertex* parNeighbor = parent->neighbors[i];
		if(AreAdjacent(parNeighbor->left, node))
			node->neighbors.push_back(parNeighbor->left);
		if(AreAdjacent(parNeighbor->right, node))
			node->neighbors.push_back(parNeighbor->right);
	}
	
	if(parent->left && (parent->left != node))
		node->neighbors.push_back(parent->left);
	
	if(parent->right && (parent->right != node))
		node->neighbors.push_back(parent->right);

	if(node->left)
		UpdateNeighbors(node->left);
	if(node->right)
		UpdateNeighbors(node->right);
	return;
}

void GetInteractionList(Vertex* node, std::vector<Vertex*>& interactionList, bool computeNeighbors=false) {
	//The list is empty for levels 1 (rootnode in this implementation) and 2 (root's children) because nodes at those levels all nodes are neighbors with each other.
	if(node->level < 3)
		return;

	//TODO (this flag is not used currently): if computeNeighbors flag is true (only in case of leaf nodes when traversing up from leaf nodes) the interaction list is all neighbors
	if(computeNeighbors)
		interactionList = node->neighbors;
	else {
		for(int i=0;i<node->neighbors.size();i++) {
			
			std::set<Vertex*> cellsConsidered;
			Vertex* neighbor = node->neighbors[i];

			if(neighbor->parent != node->parent) {
				Vertex* neighParent = neighbor->parent;
				std::pair<std::set<Vertex*>::iterator, bool> ret = cellsConsidered.insert(neighParent);
				if(ret.second) {
						if(neighParent->left && !AreAdjacent(neighParent->left,node)) {
							if((neighParent->left->radius + node->radius) <= (TAU * abs(node->center - neighParent->left->center)))
								interactionList.push_back(neighParent->left);	
							else
								node->neighbors.push_back(neighParent->left);
						}
						if(neighParent->right && !AreAdjacent(neighParent->right,node)) {
							if((neighParent->right->radius + node->radius) <= (TAU * abs(node->center - neighParent->right->center)))
								interactionList.push_back(neighParent->right);	
							else
								node->neighbors.push_back(neighParent->right);
						}
					}
				}
			}
		}
	}

void GetLeafNodes(Vertex* node) {
	if(node->isLeaf) {
		allLeaves.push_back(node);
		return;
	}
	
	if(node->left)
		GetLeafNodes(node->left);
	if(node->right)
		GetLeafNodes(node->right);
}

void TraverseDown_Recursive(Vertex* node)
{
		std::vector<Vertex*> wellSeparatedNodes;
		GetInteractionList(node, node->wellSeparatedNodes, false);
		if(node->left)
			TraverseDown_Recursive(node->left);

		if(node->right)
			TraverseDown_Recursive(node->right);
	return;
}

/* Goal: to accelerate a matrix vector multiplication Phi*q. Phi is a 'separable' matrix that can be created using column vectors x and y. 
*fun is the kernel (e.g. used to determine electrostatic potential, gravitational force etc.)
*r is the number of terms used in the Taylor series expansion. 
* org specifies the order of elements in x for determining gap.
* By default, scaling is on (scaling is used for stability) */
double* fmm1d_local_shift(int r, double *x, double *y, double * q, const double *gap, const int* org, const int fun, const int numXElems, int numYElems, const int scaling=1) {
	const double pi = 3.14159265358979323846;
	std::qsort(x, numXElems, sizeof(double), compare);  
	std::qsort(y, numYElems, sizeof(double), compare); 
	double z2 = max(x[numXElems-1], y[numYElems-1]);
	double z1 = min(x[0], y[0]);
	z2 += 0.1 * abs(z2);
	z1 -= 0.1 * abs(z1);

	Vertex* markerNode = new Vertex(); // this is a dummy node created to avoid having 'if(parent) == NULL' checks in ConstructSPatialTree.
	markerNode->level=0;

	eta0 = pow(2*pi*r, 0.5/r) / exp(1);
	Vertex* node = ConstructSpatialTree(x,y,0,numXElems,0,numYElems, (z1+z2)/2, (z2-z1)/2, 1, 1, markerNode); 
	markerNode->left=node;
	SetRootNode(node);
	

	if(node->left)
		UpdateNeighbors(node->left);

	if(rootNode->right)
		UpdateNeighbors(node->right);


	TraverseDown_Recursive(node); 
	
	return NULL;
}	
