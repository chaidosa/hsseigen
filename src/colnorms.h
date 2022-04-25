#ifndef COLNORMS_H
#define COLNORMS_H
#include<vector>
std::vector<double> colnorms(std::vector<double>& d, std::vector<double>& lam, std::vector<double>& tau, const int *org, int org_size, std::vector<double>& v, double N = 17000);
#endif
