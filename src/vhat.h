#ifndef VHAT_H
#define VHAT_H
#include <vector>
std::vector<double> vhat(std::vector<double>& d, double* lam, const int *org, int org_size,std::vector<double>& tau, std::vector<double>& w, double N = 17000);
#endif
