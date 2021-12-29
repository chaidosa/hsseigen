#include "bsxfun.h"
#include <iostream>
#include <assert.h>
using namespace std;
/*  Bsxfun similar to matlab,
*   character 'T' for times   
*/
void bsxfun(char method, double *X, std::pair<int,int>xSize, double *Y, std::pair<int,int>ySize){

    switch (method)
    {
    case 'T':
        if(xSize.first != ySize.first){
            cout<<"Error using bsxfun Non-singleton\n dimensions of the two input arrays must match each other"<<endl;
            assert(false);
        }
        for(int row = 0 ; row < xSize.first; row++){
            for(int col = 0; col<xSize.second; col++){
                X[col+row*xSize.second] = X[col+row*xSize.second]*Y[row];
            }
        }
        break;
    
    default:
        break;
    }








}