#include "bsxfun.h"
#include <iostream>
#include <assert.h>
#include <string.h>
#include <cstring>
#include <bits/stdc++.h>
#ifndef OPENBLAS
extern "C"
{
#endif
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
#ifndef OPENBLAS
}
#endif
using namespace std;
/*  Bsxfun similar to matlab,
*   character 'T' for times
*   character 'P' for plus
*   character 'M' for minus   
*/


void bsxfun(char method, double **tempX, std::pair<int,int>xSize, double *Y, std::pair<int,int>ySize)
{
    double *X = *tempX;
    switch (method)
    {
    //@times
    case 'T':
        if(xSize.first != ySize.first){
            std::cout<<"Error using bsxfun (Times) Non-singleton\n dimensions of the two input arrays must match each other\n";
            assert(false);
        }
        for(unsigned int row = 0 ; row < xSize.first; row++){
            for(unsigned int col = 0; col<xSize.second; col++){
                X[col+row*xSize.second] = X[col+row*xSize.second]*Y[row];
            }
        }
        break;
    //@plus
    case 'P':
        if(xSize.first != ySize.first){
            std::cout<<"Error using bsxfun (Plus) Non-singleton\n dimensions of the two input arrays must match each other\n";
            assert(false);
        }
        for(unsigned int row = 0 ; row < xSize.first; row++){
            for(unsigned int col = 0; col<xSize.second; col++){
                X[col+row*xSize.second] = X[col+row*xSize.second]+Y[row];
            }
        }
        break;
    //@plus row-wise: e.g. adding A+B, where A is a matrix and B is a vector would yield adding the first element of B to all the elements in the first row of A, adding second element of B to all the elements of the second row of A, and so on. So, num rows of A and B must match.
    case 'p':
        if(xSize.first != ySize.first){
            std::cout<<"Error using bsxfun (p column wise) Non-singleton\n dimensions of the two input arrays must match each other\n";
            assert(false);
        }
        for(unsigned int row = 0 ; row < xSize.first; row++){
            for(unsigned int col = 0; col<xSize.second; col++){
                X[col+row*xSize.second] += Y[row];
            }
        }
        break;

    //@minus
    case 'M':
        if(xSize.first != ySize.first){
            std::cout<<"Error using bsxfun (Minus) Non-singleton\n dimensions of the two input arrays must match each other\n";
            assert(false);
        }
        for(unsigned int row = 0 ; row < xSize.first; row++){
            for(int col = 0; col<xSize.second; col++){
                X[col+row*xSize.second] = X[col+row*xSize.second]-Y[row];
            }
        }
        break;
    //@minus column-wise: subtracting A-B, where A is a matrix and B is a vector. Num cols of A and B must match.
    case 'm':
		if(xSize.second != ySize.second){
		    std::cout<<"Error using bsxfun (Minus columnwise) Non-singleton\n dimensions of the two input arrays must match each other\n";
		    assert(false);
		}
		for(unsigned int row = 0 ; row < xSize.first; row++){
		    for(int col = 0; col<xSize.second; col++){
			X[col+row*xSize.second] = X[col+row*xSize.second]-Y[col];
		    }
		}
        break;

    default:
        std::cout<<"Enter valid Parameters";
        assert(false);
        break;

    
    }
}

//used for arranging vector
void arrange_elements(double *Arr,std::pair<int,int>arrSize,double *Indices, std::pair<int,int>indSize,double *result){
    
    if(arrSize.first != indSize.first){
        std::cout<<"Enter valid arguments";
        assert(false);
    }

    for(unsigned int i = 0;i < arrSize.first; i++){
        result[i] = Arr[(int)Indices[i]];
    }
}

//used for arranging  arrays
void arrange_elements2(double *Xc, std::pair<int,int>XSize,int *Tt, std::pair<int,int>tSize){
    if(XSize.first == 1)
        return;
    assert(XSize.first == tSize.first);
     LAPACKE_dlapmr(LAPACK_ROW_MAJOR,true, XSize.first, XSize.second, Xc, XSize.second, Tt);
}


double vec_norm(std::vector<double> &v){
    double result = 0;
    for(unsigned int i = 0; i <(int)v.size(); ++i){
        result +=v[i]*v[i];
    }

return sqrt(result);    
}

std::vector<double> diff_vec(std::vector<double> V){
    std::vector<double> result;
    for(unsigned int i = 0; i <(int)V.size()-1; ++i){
        double temp = V[i+1]-V[i];
        result.push_back(temp);
    }
return result;
}
