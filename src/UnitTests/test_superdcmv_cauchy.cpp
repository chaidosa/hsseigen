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
    }
}

template<typename T>
void loadXX(int Xsize, ifstream& fp, T* X){
    setprecision(32);
    for (int i = 0; i < Xsize; i++)
    {   
        T temp;
        fp >>  temp;
        X[i] = temp-1;
    }
}

int main(){
    ifstream input("/workspaces/hsseigen/src/UnitTests/SuperdcmvCauchydump.txt");
    pair<int,int> Qsize;
    input >> Qsize.first; input >> Qsize.second;

    pair<int, int> Qcsize;
    input >> Qcsize.first; input >> Qcsize.second;

    pair<int, int> Qcsizes[6];
    input >> Qcsizes[0].first; input >> Qcsizes[0].second;
    input >> Qcsizes[1].first; input >> Qcsizes[1].second;
    input >> Qcsizes[2].first; input >> Qcsizes[2].second;
    input >> Qcsizes[3].first; input >> Qcsizes[3].second;
    input >> Qcsizes[4].first; input >> Qcsizes[4].second;
    input >> Qcsizes[5].first; input >> Qcsizes[5].second;

    pair<int, int> Jsize;
    input >> Jsize.first; input >> Jsize.second; 
    
    pair<int, int> Gsize;
    input >> Gsize.first; input >> Gsize.second;

    pair<int, int> Isize;
    input >> Isize.first; input >> Isize.second;

    pair<int, int> v2csize;
    input >> v2csize.first; input >> v2csize.second;

    pair<int, int> Tsize;
    input >> Tsize.first; input >> Tsize.second;

    pair<int, int> Xsize;
    input >> Xsize.first; input >> Xsize.second;

    double *QC[5]; int *org;
    for (int i = 0; i < 5; i++)
    {
        QC[i] = new double[Qcsizes[i].first]; loadX(Qcsizes[i].first, input, QC[i]);
    }
    org = new int[Qcsizes[5].first]; loadXX(Qcsizes[5].first, input, org);
    
    double *J = new double[Jsize.first]; loadXX(Jsize.first, input, J);
    double *G = new double[Gsize.first]; loadX(Gsize.first, input, G); G[0] -= 1; G[1] -=1;
    double *I = new double[Isize.first]; loadXX(Isize.first, input, I);
    double *v2c = new double[v2csize.first]; loadX(v2csize.first, input, v2c);
    int *T = new int[Tsize.first]; loadXX(Tsize.first, input, T);
    double *X = new double[Xsize.first]; loadX(Xsize.first, input, X);

    int n, n1, n2, n3;
    input >> n;
    input >> n1; input >> n2; input >> n3;

    double *Xf = new double[Xsize.first]; loadX(Xsize.first, input, Xf);
    
    nonleaf *Qq = new nonleaf();

    for (int i = 0; i < 5; i++)
    {
         Qq->QC[i] = QC[i];
    }
    
    for (int i = 0; i < 6; i++)
    {
        Qq->qcSizes[i] = Qcsizes[i];    
    }
    
    Qq->Org = org;
    Qq->J = J;
    Qq->G = G;
    Qq->I = I;
    Qq->v2c = v2c;
    Qq->T = T;
    Qq->n = n; Qq->n1 = n1; Qq->n2 = n2; Qq->n3 = n3;

    Qq->JSize = Jsize;
    Qq->GSize = Gsize;
    Qq->ISize = Isize;
    Qq->v2cSize = v2csize;
    Qq->TSize = Tsize;

    superdcmv_cauchy(Qq, Qsize, X, Xsize, 0, 1024);

    double sqrSum = 0;

    for (int i = 0; i < 31; i++)
    {
        // cout << (abs(Yexpected[i] - Y[i]) < 1e-10) << "\n";
    }
    

    for (int i=0; i<Xsize.first; i++){
        cout << "X, Xexpected " << X[i] << " " << Xf[i] << "\n";
        sqrSum += (X[i]-Xf[i])*(X[i]-Xf[i]);
    }

    cout << "Rms error is: " << sqrt(sqrSum/Xsize.first) << "\n";

    input.close();
}