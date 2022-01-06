#ifndef BSXFUN_H
#define BSXFUN_H
#include <string>
#include <iostream>
#include <assert.h>
#include <cstring>
void bsxfun(char method, double **X, std::pair<int,int>xSize, double *Y, std::pair<int,int>ySize);
void arrange_elements(double *Arr,std::pair<int,int>arrSize,double *Indices, std::pair<int,int>indSize,double *result);
//used for arrange  arrays
void arrange_elements2(double *X, std::pair<int,int>XSize,double *T, std::pair<int,int>tSize,double *result);
#endif