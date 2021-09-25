#ifndef DIVIDE_H
#define DIVIDE_H
#include"mat2hsssym.h"
typedef struct DivideOutParams
{
	std::vector<int> l;
	std::vector<std::pair<int, int> > lvec, dSize;
	double* Z, **D;
	std::pair<int, int> zSize;
}DivideOutParams;

DivideOutParams* Divide(HSSMat* hssMat, BinTree* tr);


#endif
