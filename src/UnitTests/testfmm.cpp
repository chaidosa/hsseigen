#include<iostream>
#include<fstream>
#include<iomanip>
#include<cassert>
#include"../fmm1d_local_shift_2.h"
#include"../fmm_types.h"

//the input file trifmminput_ex4_4k_fun1.txt contains x,y,number of terms, gap, q, and org in that order for the input ex4_4k.txt used in the first call to trifmm1d_local_shift with function 1. in rootfinder.m of MATLAB code.
//the input file trifmminput_ex4_4k_fun2.txt contains x,y,number of terms, gap, q, and org in that order for the input ex4_4k.txt used in the third call to trifmm1d_local_shift with function 2. in rootfinder.m of MATLAB code.
void PrintTree(const Vertex* node);

using namespace std;
int main(int argc, char* argv[]){
	if(argc!=2){
		printf("usage:./fmm testfilename\n");
		return 0;
	}
	string header;
	int numXElems=0, numYElems=0, numQElems=0, numR=0, numGapElems=0, numOrgElems=0; 
	ifstream ifs(argv[1],ifstream::in);
	if(!ifs) {
		cerr<<"Error: unable to open the file fmminput_out_32.txt. exiting."<<endl;
		exit(1);
	}

	getline(ifs,header);
	ifs>>numR>>numXElems>>numYElems>>numGapElems>>numQElems>>numOrgElems;
	cout<<header<<endl;
	cout<<"header: numXElems("<<numXElems<<"), numYElems("<<numYElems<<"), numR("<<numR<<"), numGapElems("<<numGapElems<<"), numQElems("<<numQElems<<"), numOrgElems("<<numOrgElems<<endl;

	assert(numXElems == numOrgElems);
	double* x= new double[numXElems];
	double* y= new double[numYElems];
	double* q= new double[numQElems];
	double* gap= new double[numGapElems];
	int* org= new int[numOrgElems];
	int* r=new int[numR];



	for(int i=0;i<numR;i++)
		ifs>>r[i];
	for(int i=0;i<numXElems;i++)
		ifs>>x[i];

	for(int i=0;i<numYElems;i++)
		ifs>>y[i];
	
	for(int i=0;i<numGapElems;i++)
		ifs>>gap[i];
	
	for(int i=0;i<numQElems;i++)
		ifs>>q[i];

	for(int i=0;i<numXElems;i++)
		ifs>>org[i];

	cout<<setprecision(12)<<"x["<<numXElems-1<<"]="<<x[numXElems-1]<<endl;
	cout<<setprecision(12)<<"y["<<numYElems-1<<"]="<<y[numYElems-1]<<endl;
	cout<<"r["<<numR-1<<"]="<<r[numR-1]<<endl;
	cout<<setprecision(12)<<"gap["<<numGapElems-1<<"]="<<gap[numGapElems-1]<<endl;
	cout<<setprecision(12)<<"q["<<numQElems-1<<"]="<<q[numQElems-1]<<endl;
	cout<<setprecision(12)<<"org["<<numOrgElems-1<<"]="<<org[numOrgElems-1]<<endl;
	ifs.close();

	//double* z1=fmm1d_local_shift_2(r[0],x,y,q,gap,org,3,numXElems,numYElems);
	//double* z2=fmm1d_local_shift(r[0],x,y,q,gap,org,3,numXElems,numYElems);
	double* z3=trifmm1d_local_shift(r[0],x,y,q,gap,org,2,numXElems,numYElems);

	delete [] x;
	delete [] y;
	delete [] q;
	delete [] gap;
	delete [] org;
	delete [] r;
}

void PrintTree(const Vertex* node) {
	if(node == NULL)
		return;
	printf("label:%d\t(c, r):(%lf %lf)\tlvl:%d\n",node->label, node->center, node->radius, node->level);
	printf("\tnbrs:");
	for(int i=0;i<node->nbrs.size();i++)
		printf("%d ",node->nbrs[i]);
	printf("\n");
	printf("\til:");
	for(int i=0;i<node->il.size();i++)
		printf("%d ",node->il[i]);
	printf("\n");
	printf("\t(numX, numY):(%d, %d)\t(xL,xR,yL,yR):(%d %d %d %d)\n",node->xRight-node->xLeft,node->yRight-node->yLeft, node->xLeft, node->xRight, node->yLeft, node->yRight);
	PrintTree(node->left);
	PrintTree(node->right);
}

