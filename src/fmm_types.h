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
	int resSize;
#endif
	short int level;
	bool isLeaf;
	int xLeft, xRight, yLeft, yRight; //indices of points within the box
	double center, radius;//, eta;
	double* res; //used to hold intermediate results v
	double* ul, *uu, *zl, *zu; //used to hold intermediate results u, and z. In case of triangular fmm zl and zu are produced. z is stored in zl and u in ul in case of normal fmm with shift.
	bool computed, computedUpper;
	Vertex():parent(0), level(0), isLeaf(false), left(NULL), right(NULL){
		xLeft=-1;xRight=-1;yLeft=-1;yRight=-1;res=NULL;computed=false;computedUpper=false; ul=NULL;uu=NULL; zl=NULL; zu=NULL;
#ifdef DEBUG
		label=0;
		resSize=0;
#endif
	}
	~Vertex(){
		if(res){
			delete [] res;
			res=NULL;
		}
		if(ul) {
			delete [] ul;
			ul=NULL;
		}
		if(uu) {
			delete [] uu;
			uu=NULL;
		}
		if(zl) {
			delete [] zl;
			zl=NULL;
		}
		if(zu) {
			delete [] zu;
			zu=NULL;
		}

	}
}Vertex;



#ifdef DEBUG
void AssignLabels(Vertex* node);
#endif

#endif

