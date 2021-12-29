#ifndef BSXFUN_H
#define BSXFUN_H
#include <string>
#include <istream>
void bsxfun(char method, double *X, std::pair<int,int>xSize, double *Y, std::pair<int,int>ySize);
void arrange_elements(double *Arr,std::pair<int,int>arrSize,double *Indices, std::pair<int,int>indSize,double *result);
#endif