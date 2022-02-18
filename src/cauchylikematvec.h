#ifndef CAUCHYLIKEMATVEC_H
#define CAUCHYLIKEMATVEC_H
#include <utility>
void cauchylikematvec(double ***Qcd, std::pair<int,int>*qcSizes,double **Xx, std::pair<int,int>xSize,int ifTrans, double N=1024);
#endif