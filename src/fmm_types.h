#ifndef __FMM_TYPES_H
#define __FMM_TYPES_H

#define MAX_POINTS_IN_CELL 512
#define MAX_LEVELS 64
#include<vector>

typedef struct Vertex{
	Vertex* left, *right;
	Vertex* parent;
	std::vector<Vertex*> neighbors;
	std::vector<Vertex*> wellSeparatedNodes;
	int label;
	short int level;
	bool isLeaf;
	int xLeft, xRight, yLeft, yRight; //indices of points within the box
	double center, radius, eta;
	Vertex():parent(0), level(0), label(0), isLeaf(false), left(NULL), right(NULL){xLeft=-1;xRight=-1;yLeft=-1;yRight=-1;}
}Vertex;


const Vertex* GetRootNode();
void SetRootNode(Vertex* node);

#endif

