#ifndef __FMM_TYPES_H
#define __FMM_TYPES_H

#define MAX_POINTS_IN_CELL 512
#define MAX_LEVELS 64
#include<vector>
#include<float.h>
#include<limits.h>

typedef struct Vertex{
	Vertex* left, *right;
	Vertex* parent;
	std::vector<Vertex*> neighbors;
	long int label;
	short int level;
	bool isLeaf;
	int xLeft, xRight, yLeft, yRight; //indices of points within the box
	Vertex():parent(0), level(0), label(0), isLeaf(false), left(NULL), right(NULL){}
}Vertex;


#endif

