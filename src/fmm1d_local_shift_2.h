#ifndef FMM1D_LOCAL_SHIFT
#define FMM1D_LOCAL_SHIFT
double* fmm1d_local_shift_2(int r, double *x, double *y, double * q, const double *gap, const int* org, const int fun, const int numXElems, int numYElems, const int scaling=1);
double* fmm1d_local_shift(int r, double *x, double *y, double * q, const double *gap, const int* org, const int fun, const int numXElems, int numYElems, const int scaling=1);
double* trifmm1d_local_shift(int r, double *x, double *y, double * q, const double *p_gap, const int* p_org, const int fun, const int numXElems, int numYElems, const int scaling=1);
#endif
