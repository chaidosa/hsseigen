#include "bsxfun.h"
#include <iostream>
#include <assert.h>
#include <string.h>
#include <cstring>
#include <bits/stdc++.h>
using namespace std;
/*  Bsxfun similar to matlab,
*   character 'T' for times
*   character 'P' for plus
*   character 'M' for minus   
*/
void bsxfun(char method, double **tempX, std::pair<int,int>xSize, double *Y, std::pair<int,int>ySize){
    double *X = *tempX;
    switch (method)
    {
    //@times
    case 'T':
        if(xSize.first != ySize.first){
            std::cout<<"Error using bsxfun Non-singleton\n dimensions of the two input arrays must match each other\n";
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
            std::cout<<"Error using bsxfun Non-singleton\n dimensions of the two input arrays must match each other\n";
            assert(false);
        }
        for(unsigned int row = 0 ; row < xSize.first; row++){
            for(unsigned int col = 0; col<xSize.second; col++){
                X[col+row*xSize.second] = X[col+row*xSize.second]+Y[row];
            }
        }
        break;
    //@minus
    case 'M':
        if(xSize.first != ySize.first){
            std::cout<<"Error using bsxfun Non-singleton\n dimensions of the two input arrays must match each other\n";
            assert(false);
        }
        for(unsigned int row = 0 ; row < xSize.first; row++){
            for(int col = 0; col<xSize.second; col++){
                X[col+row*xSize.second] = X[col+row*xSize.second]-Y[row];
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
void arrange_elements2(double *X, std::pair<int,int>XSize,double *T, std::pair<int,int>tSize,double *result){
    double *temp = new double[(tSize.second * XSize.second)];

    for(unsigned int row = 0; row < tSize.second; row++){
        memcpy(temp+row*(XSize.second),X+((int)T[row])*(XSize.second),sizeof(double)*(XSize.second));
    }
    result = temp;
    temp = NULL;
}


double vec_norm(std::vector<double> v){
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