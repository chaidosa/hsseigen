#ifndef SUPERDCMV_CAUCHY_H
#define SUPERDCMV_CAUCHY_H
#include <iostream>
#include <utility>
void superdcmv_cauchy(double **Q,std::pair<int, int>*qSize, double *X,std::pair<int, int>xSize,int ifTrans,double N = 1024);


#endif