#ifndef __FMM_TYPES_H
#define __FMM_TYPES_H

#define MAX_POINTS_IN_CELL 512 //max. number of points in leaf-node
#define MAX_LEVELS 64
#include<vector>
#ifdef DEBUG
#include<queue>
#endif

typedef struct Vertex{
	Vertex* left, *right;
	Vertex* parent;
	std::vector<Vertex*> neighbors;
	std::vector<Vertex*> wellSeparatedNodes;
#ifdef DEBUG
	int label;
	std::vector<int> nbrs, il;
	int resuSize;
#endif
	short int level;
	bool isLeaf;
	int xLeft, xRight, yLeft, yRight; //indices of points within the box
	double center, radius;//, eta;
	double* res, *resu; //used to hold intermediate results v, u, and z. resu is used to hold uppertriangular result of u and z
	bool computed, computedUpper;
	Vertex():parent(0), level(0), isLeaf(false), left(NULL), right(NULL){
		xLeft=-1;xRight=-1;yLeft=-1;yRight=-1;res=NULL;resu=NULL;computed=false;computedUpper=false;
#ifdef DEBUG
		label=0;
		resuSize=0;
#endif
	}
}Vertex;


#ifdef DEBUG
void AssignLabels(Vertex* node);
#endif
const Vertex* GetRootNode();
void SetRootNode(Vertex* node);

#endif

