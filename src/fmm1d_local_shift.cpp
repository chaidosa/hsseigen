#include<cstdlib>
#include<cassert>
#include"fmm_types.h"
#define PARTITION_THRESHOLD 512 //number of points in leaf-node
#define SEPARATION_RATIO 0.6
/* Goal: to a matrix vector multiplication Phi*q. Phi is a 'separable' matrix that can be created using column vectors x and y. 
*fun is the kernel (e.g. used to determine electrostatic potential, gravitational force etc.)
*r is the number of terms used in the Taylor series expansion. 
* org specifies the order of elements in x for determining gap.
* By default, scaling is on (scaling is used for stability) */
using namespace std;
typedef pair<double*, int*> TreeNodeData;

int compare(const void* x, const void* y) {
	const double arg1 = *static_cast<const double*>(x);
	const double arg2 = *static_cast<const double*>(y);
	if (arg1 < arg2) return -1;
	if (arg1 > arg2) return 1;
	return 0;
}

Vertex* ConstructSpatialTree(double *x, double *y, int xLb, int xUb, int yLb, int yUb, double center, double rad, int el, int er, Vertex* parent) {

	double c = center - rad/2; 
	double d = rad/2;
	double mid;
	bool flag=false;
	int i, j, k, l, countXLeft=0, countXRight=0, countYLeft=0, countYRight=0;
	Vertex* leftChild=NULL, *rightChild=NULL, *node=NULL;

	if((xLb > xUb) || (yLb > yUb))
		return NULL;

	int sizex = xUb-xLb, sizey = yUb-yLb;
	if ((sizex <= MAX_POINTS_IN_CELL) || (sizey <= MAX_POINTS_IN_CELL)){
		node = new Vertex();
		node->xLeft = xLb;
		node->xRight = xUb;
		node->yLeft = yLb;
		node->yRight = yUb;
		node->isLeaf = true;
		node->level = parent->level + 1;
		return node;
	}

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
		//x and y intervals are empty for the bisected point. Adjust center and rad.
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

	return node;
}

double* fmm1d_local_shift(int r, double *x, double *y, double * q, const double *gap, const int* org, const int fun, const int numXElems, int numYElems, const int scaling=1) {
	
	std::qsort(x, numXElems, sizeof(double), compare);  
	std::qsort(y, numYElems, sizeof(double), compare); 
	double z2 = max(x[numXElems-1], y[numYElems-1]);
	double z1 = min(x[0], y[0]);
	z2 += 0.1 * abs(z2);
	z1 -= 0.1 * abs(z1);

	Vertex* markerNode = new Vertex(); // this is a dummy node created to avoid having 'if(parent) == NULL' checks in ConstructSPatialTree.
	markerNode->level=0;

	Vertex* rootNode = ConstructSpatialTree(x,y,0,numXElems,0,numYElems, (z1+z2)/2, (z2-z1)/2, 1, 1, markerNode); 
	return NULL;
}	
