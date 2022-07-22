/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% compute vhat by Lowner's formula %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include<bits/stdc++.h>
#include<fstream>
#include<iomanip>
#include <iostream>
#include<fstream>
#include<iomanip>
#include "fmm1d_local_shift_2.h"
#include "bsxfun.h"
#include "vhat.h"

/*template<typename T>
void PrintArray(T *A, int row, int col, const char* filename="debugoutput.txt")
{
#if 1
    ofstream txtOut;
    txtOut.open(filename, std::ofstream::out | std::ofstream::app);
    for(unsigned int i = 0; i <(unsigned)row; i++)
    {
        for(int j=0; j < col; j++)
        {
            txtOut<<setprecision(12)<<A[j+i*col]<<"\n";
        }
    }
    txtOut.close();
#endif
}*/


double * vhat(std::vector<double>& d, double* lam, const int *org, int org_size, std::vector<double>& tau, std::vector<double>& w, double N){
    /*ofstream txtOut;
    txtOut.open("vhat_input.txt", std::ofstream::out | std::ofstream::app);
    txtOut<<d.size()<<"\n"<<org_size<<"\n"<<org_size<<"\n"<<tau.size()<<"\n"<<w.size()<<"\n"<<N<<"\n";
    for(int i=0;i<d.size();i++)
	txtOut<<setprecision(12)<<d[i]<<"\n";
    for(int i=0;i<org_size;i++)
	txtOut<<setprecision(12)<<lam[i]<<"\n";
    for(int i=0;i<org_size;i++)
	txtOut<<setprecision(12)<<org[i]<<"\n";
    for(int i=0;i<tau.size();i++)
	txtOut<<setprecision(12)<<tau[i]<<"\n";
    for(int i=0;i<w.size();i++)
	txtOut<<setprecision(12)<<w[i]<<"\n";
    txtOut.close();*/

    int n = org_size;
    int r = 50;   
    //std::vector<double>v(n);
    double *v = new double[n];
  //  if(n < N) 
    {
        int dRows = d.size();
        int dCols = org_size;
        for(int row = 0; row < dRows; row++){
            v[row] = 0;
            for(int col = 0; col < dCols; col++){
                double temp = d[row] - d[org[col]];
                temp = temp - tau[col];
                //Dlam[col + row*dCols]
                //D[col + row*dRows] = d[row] - d[col];
                double temp2 = d[row] - d[col];
                if(row == col)
                    temp2 = 1;

                v[row] += std::log(std::abs(temp)) -std::log(std::abs(temp2));

            }

            double temp3 = v[row]/2;
            v[row] = std::exp(temp3);
        }
        for(int itr = 0; itr < w.size(); itr++){
            if(w[itr] < 0)
                v[itr] = -v[itr];
        }              
    }
#if 0
    else
    {
	assert(d.size() == n);
	assert(w.size() == n);
        
	double* ptrD=d.data();
	double* gap=tau.data();
	std::vector<double> e1(n,1);
	std::vector<double> zeros(n,0);
	
	/*char fname[64];
	sprintf(fname,"vhat_fmminput_call1.txt");
	std::ofstream txtOut;
	txtOut.open(fname, std::ofstream::out | std::ofstream::app);
	txtOut<<"50\n";// dummy. to satify matlab load("") command.
	txtOut<<"50\n";// dummy. to satify matlab load("") command.
	txtOut<<r<<"\n";
	txtOut.close();
	PrintArray<double>(ptrD,1,d.size(),fname);
	PrintArray<double>(lam,1,n,fname);
	double* tempe1=e1.data();
	PrintArray<double>(tempe1,1,n,fname);
	PrintArray<double>(gap,1,d.size(),fname);
	PrintArray<const int>(org,1,d.size(),fname);
	txtOut.close();*/

	double* v0 = fmm1d_local_shift_2(r, ptrD, lam, e1.data(), gap, org, 3, d.size(), n);

	std::vector<int> tmpOrg(d.size(),0);
	for(int i=0;i<d.size();i++)
		tmpOrg[i]=i;

	/*sprintf(fname,"vhat_fmmoutput_call1_v0.txt");
	PrintArray<double>(v0,1,d.size(),fname);
	
	sprintf(fname,"vhat_fmminput_call2.txt");
	txtOut.open(fname, std::ofstream::out | std::ofstream::app);
	txtOut<<"50\n";// dummy. to satify matlab load("") command.
	txtOut<<"50\n";// dummy. to satify matlab load("") command.
	txtOut<<r<<"\n";
	txtOut.close();
	PrintArray<double>(ptrD,1,d.size(),fname);
	PrintArray<double>(ptrD,1,d.size(),fname);
	tempe1=e1.data();
	PrintArray<double>(tempe1,1,n,fname);
	tempe1=zeros.data();
	PrintArray<double>(tempe1,1,d.size(),fname);
	const int* tempe2=tmpOrg.data();
	PrintArray<const int>(tempe2,1,d.size(),fname);
	txtOut.close();*/

        double* vd = fmm1d_local_shift(r, ptrD, ptrD, e1.data(), zeros.data(), tmpOrg.data(), 3, d.size(), d.size());
	/*sprintf(fname,"vhat_fmmoutput_call2_vd.txt");
	PrintArray<double>(vd,1,d.size(),fname);*/
	for(int i=0;i<n;i++)
		v[i]=exp(0.5 *(v0[i]-vd[i]));
 
	//Vec:
	for(int i=0;i<n;i++)
		if(w[i] <0)
			v[i] = -1 * v[i];
	/*sprintf(fname,"vhat_fmmoutput_call2_v.txt");
	PrintArray<double>(v,1,n,fname);*/

    	delete [] v0;
	delete [] vd;
    }
#endif
    
    /*txtOut.open("vhat_output.txt", std::ofstream::out | std::ofstream::app);
    txtOut<<setprecision(12)<<d.size()<<"\n";
    for(int i=0;i<d.size();i++)
	txtOut<<setprecision(12)<<v[i]<<"\n";
    txtOut.close();*/

    return v;
}
