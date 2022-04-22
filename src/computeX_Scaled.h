#ifndef COMPUTEX_SCALED_H
#define COMPUTEX_SCALED_H
void ComputeU_Scaled(double** UTrans, const double* x, int numXElems, int r, double eta, double a, double dx, const int scaling);
void ComputeB_Scaled(double** B, int r, double eta0, double a, double b, double dx, double dy, int fun, const int scaling);
void ComputeT_Scaled(double** T, int r, double eta0, double a, double b, double dx, double dy, const int scaling);
#endif
