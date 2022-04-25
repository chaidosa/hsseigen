#ifndef EIGENMATRIX_H
#define EIGENMATRIX_H
#include<utility>
#include<iostream>
class nonleaf{
        public:
        double *QC[5];
        std::pair<int,int>qcSizes[6];        
        int *Org;
        
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

        ~nonleaf(){
                delete[] QC[0];
                delete[] QC[1];
                delete[] QC[2];
                delete[] QC[3];
                delete[] QC[4];
                delete[] Org;
                delete[] J;
                delete[] G;
                delete[] I;
                delete[] v2c;
                delete[] T;   
                std::cout<<"deleted non Leaf\n";                             

        }
}; 
class EIG_MAT{
    public:
    double  *Q0_leaf;    
    nonleaf **Q0_nonleaf;
    int n_non_leaf;

    ~EIG_MAT(){
        delete[] Q0_leaf;        
    }    
};
	
#endif
