#ifndef COLNORMS_H
#define COLNORMS_H
#include<vector>
double* colnorms(std::vector<double>& d, double* lam, std::vector<double>& tau, const int *org, int org_size, std::vector<double>& v, double N = 1024);
#endif
