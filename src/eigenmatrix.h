#ifndef EIGENMATRIX_H
#define EIGENMATRIX_H
#include<utility>
class nonleaf{
        public:
        double *QC[6];
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

/* Nikhil:
 * from secular.m, Q, the output of secular.m is always a 7x1 cell. Further, Q{1}, Q3, is always a 6x1 cell. So, how about the following organization? you can represent a vector using the same class. 
 * 
   class EIGMatrix{
   public: 
   	double *Q3[6]; //Q3[0], v_hat, can be represented by a double pointer. Similarly, Q3[1], s3, can be represented by a double pointer and so on..
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
     // We will have N = number of nodes in the tree, number of EIGMatrix objects. For leaf nodes, Q3[0] - Q3[5], J, G, I, v2c will be NULL and n2, n3 will always be 0. For non-leaf nodes, you need to allocate memory.
     // If wasting space for leaf nodes to hold J, G, I, v2c bothers you, you could define a parent class that don't have these fields (and then subclass it for non-leaf nodes).  
  */

	
#endif
