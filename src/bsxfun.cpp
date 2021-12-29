#include "bsxfun.h"
#include <iostream>
#include <assert.h>
using namespace std;
/*  Bsxfun similar to matlab,
*   character 'T' for times
*   character 'P' for plus
*   character 'M' for minus   
*/
void bsxfun(char method, double *X, std::pair<int,int>xSize, double *Y, std::pair<int,int>ySize){

    switch (method)
    {
    //@times
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
    //@plus
    case 'P':
        if(xSize.first != ySize.first){
            cout<<"Error using bsxfun Non-singleton\n dimensions of the two input arrays must match each other"<<endl;
            assert(false);
        }
        for(int row = 0 ; row < xSize.first; row++){
            for(int col = 0; col<xSize.second; col++){
                X[col+row*xSize.second] = X[col+row*xSize.second]+Y[row];
            }
        }
        break;
    //@minus
    case 'M':
        if(xSize.first != ySize.first){
            cout<<"Error using bsxfun Non-singleton\n dimensions of the two input arrays must match each other"<<endl;
            assert(false);
        }
        for(int row = 0 ; row < xSize.first; row++){
            for(int col = 0; col<xSize.second; col++){
                X[col+row*xSize.second] = X[col+row*xSize.second]-Y[row];
            }
        }
        break;
    default:
        cout<<"Enter valid Parameters";
        assert(false);
        break;
    }
}
void arrange_elements(double *Arr,std::pair<int,int>arrSize,double *Indices, std::pair<int,int>indSize,double *result){
    if(arrSize.first != indSize.first){
        cout<<"Enter valid arguments";
        assert(false);
    }

    for(int i = 0;i < arrSize.first; i++){
        result[i] = Arr[(int)Indices[i]];
    }
}