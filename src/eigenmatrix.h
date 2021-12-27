#ifndef EIGENMATRIX_H
#define EIGENMATRIX_H
#include<utility>
class nonleaf{
        public:
        double **QC;
        std::pair<int,int>*qcSizes;
        double *J;
        std::pair<int,int>JSize;
        double *G;
        std::pair<int,int>GSize;
        double *I;
        std::pair<int,int>ISize;
        double *v2c;
        std::pair<int,int>v2cSize;
        double *T;
        std::pair<int,int>TSize;
        int n,n1,n2,n3;
}; 
class EIG_MAT{
    public:
    double  *Q0_leaf;    
    nonleaf **Q0_nonleaf;
    std::pair<int,int>*q0_nonleaf_Sizes;
};
#endif