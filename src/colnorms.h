#ifndef COLNORMS_H
#define COLNORMS_H
#include<vector>
double* colnorms(std::vector<double>& d, double* lam, std::vector<double>& tau, const int *org, int org_size, double *v, double N = 17000);
#endif
