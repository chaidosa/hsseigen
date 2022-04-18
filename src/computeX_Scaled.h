#ifndef COMPUTEX_SCALED_H
#define COMPUTEX_SCALED_H
void ComputeU_Scaled(double* x, double eta, double a, double dx, const int scaling);
void ComputeB_Scaled(double** UTrans, double* x, int numXElems, int r, double etax, double etay, double a, double b, double dx, double dy, int fun, const int scaling);
void ComputeT_Scaled(double** UTrans, double* x, int numXElems, int r, double eta0, double a, double b, double dx, double dy, const int scaling);
#endif
