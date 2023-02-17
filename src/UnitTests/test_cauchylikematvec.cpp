#include "../cauchylikematvec.h"
#include "../bsxfun.h"
#include <string.h>
#include "../fmm1d_local_shift_2.h"
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

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

template<typename T>
void loadX(int Xsize, ifstream& fp, T* X){
    for (int i = 0; i < Xsize; i++)
    {
        fp >> X[i];
    }
}


int main(int argc, char**argv){
    ifstream inputfile("cauchylikematvec.txt");

    pair<int,int>* qcSize = new pair<int,int>[6];
    pair<int,int>  xsize;

    inputfile >> qcSize[5].first;

    inputfile >> qcSize[2].first;

    inputfile >> xsize.second; xsize.first = 31;

    qcSize[0].first = 31; qcSize[0].second = 1;
    qcSize[1].first = 31; qcSize[1].second = 1;
    qcSize[2].second = 1;
    qcSize[3].first = 31; qcSize[3].second = 1;
    qcSize[4].first = 31; qcSize[4].second = 1;
    qcSize[5].second = 1;


    int vSize = qcSize[0].first;
    double* V = new double[vSize];
    loadX(vSize, inputfile, V);

    
    int sSize = qcSize[1].first;
    double* S = new double[sSize];
    loadX(sSize, inputfile, S);

    
    int dSize = qcSize[2].first;
    double* D = new double[dSize];
    loadX(dSize, inputfile, D);

    
    int lamSize = qcSize[3].first;
    double* Lam = new double[lamSize];
    loadX(lamSize, inputfile, Lam);

    
    int tauSize = qcSize[4].first;
    double* Tau = new double[tauSize];
    loadX(tauSize, inputfile, Tau);

    
    int orgSize = qcSize[5].first;
    int* Org = new int[orgSize];
    loadX(orgSize, inputfile, Org);

    
    int xSize = xsize.first;
    double* X = new double[xSize];
    loadX(xSize, inputfile, X);


    double** qc = new double*[6];
    qc[0] = V; qc[1] = S; qc[2] = D; qc[3] = Lam; qc[4] = Tau;


    double* Y = cauchylikematvec(qc, qcSize, Org, X, xsize, 1, 1024);

    double*Yexpected = new double[xsize.first];
    ifstream Ymatlab("cauchylikematvel1.txt");
    loadX(xsize.first, Ymatlab, Yexpected);

    double sqrSum = 0;

    for (int i=0; i<xsize.first; i++){
        cout << "Y, Yexpected " << Y[i] << " " << Yexpected[i] << "\n";
        sqrSum += (Y[i]-Yexpected[i])*(Y[i]-Yexpected[i]);
    }

    cout << "Rms error is: " << sqrt(sqrSum/xsize.first) << "\n";
}