#include "../superdcmv_cauchy.h"
#include <iostream>
#include <utility>
#include <string.h>
#include "../eigenmatrix.h"
#include "../bsxfun.h"
#include "../cauchylikematvec.h"
#include <fstream>
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

#include <iomanip>

template<typename T>
void loadX(int Xsize, ifstream& fp, T* X){
    setprecision(32);
    for (int i = 0; i < Xsize; i++)
    {
        fp >>  X[i];
        fpcheck++;
    }
}

int main(){

}