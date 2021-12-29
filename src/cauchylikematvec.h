#ifndef CAUCHYLIKEMATVEC_H
#define CAUCHYLIKEMATVEC_H
#include <utility>
void cauchylikematvec(double **Qc, std::pair<int,int>*qcSizes,double *X, std::pair<int,int>xSize,int ifTrans, double N=1024);
#endif