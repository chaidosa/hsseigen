#include<iostream>
#include<fstream>
#include<iomanip>
#include"../fmm1d_local_shift.h"
using namespace std;
int main(int argc, char* argv[]){
	string header;
	int numXElems=0, numYElems=0, numQElems=0, numR=0, numGapElems=0; 
	ifstream ifs("fmminput_out_32.txt",ifstream::in);
	if(!ifs) {
		cerr<<"Error: unable to open the file fmminput_out_32.txt. exiting."<<endl;
		exit(1);
	}

	getline(ifs,header);
	ifs>>numXElems>>numYElems>>numR>>numGapElems>>numQElems;
	cout<<header<<endl;
	cout<<"header: numXElems("<<numXElems<<"), numYElems("<<numYElems<<"), numR("<<numR<<"), numGapElems("<<numGapElems<<"), numQElems("<<numQElems<<")"<<endl;

	double* x= new double[numXElems];
	double* y= new double[numYElems];
	double* q= new double[numQElems];
	double* gap= new double[numGapElems];
	int* org= new int[numXElems];
	int* r=new int[numR];


	for(int i=0;i<numXElems;i++)
		org[i]=i;

	for(int i=0;i<numXElems;i++)
		ifs>>x[i];

	for(int i=0;i<numYElems;i++)
		ifs>>y[i];
	
	for(int i=0;i<numR;i++)
		ifs>>r[i];
	
	for(int i=0;i<numGapElems;i++)
		ifs>>gap[i];
	
	for(int i=0;i<numQElems;i++)
		ifs>>q[i];

	cout<<setprecision(12)<<"x["<<numXElems-1<<"]="<<x[numXElems-1]<<endl;
	cout<<setprecision(12)<<"y["<<numYElems-1<<"]="<<y[numYElems-1]<<endl;
	cout<<"r["<<numR-1<<"]="<<r[numR-1]<<endl;
	cout<<setprecision(12)<<"gap["<<numGapElems-1<<"]="<<gap[numGapElems-1]<<endl;
	cout<<setprecision(12)<<"q["<<numQElems-1<<"]="<<q[numQElems-1]<<endl;
	ifs.close();

	fmm1d_local_shift(r[0],x,y,q,gap,org,3,numXElems,numYElems);

	delete [] x;
	delete [] y;
	delete [] q;
	delete [] gap;
	delete [] org;
	delete [] r;
}
