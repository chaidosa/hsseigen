#include<cstdlib>
#include<cstdio>
#include<cassert>
#include<math.h>
#include<cstring>
#include"computeX_Scaled.h"
using namespace std;

//function computes local expansion.
//x is the vector of elements in the cell, r is the number of terms in taylor expansion, eta is a per-node term computed as ((2*pi*r)^(0.5/r))/exp(1), a is the center of the cell, dx is the diameter of the cell.
void ComputeU_Scaled(double** UTrans, const double* x, int numXElems, int r, double eta0, double a, double dx, const int scaling) {

	*UTrans = new double[numXElems*r];
	double* temp=*UTrans;
	double eta = eta0*2/dx;
		if(scaling == 1) {
		
		//1st term in the multipole expansion ia 1.
		for(int i=0;i<numXElems;i++) temp[i]=1;
		
		//2nd term in the multipole expansion
		for(int i=0;i<numXElems;i++)
			temp[i+numXElems]=eta * (x[i] - a); 

		for(int k=1;k<r-1;k++)
			for(int i=0;i<numXElems;i++)
				temp[numXElems*(k+1)+i]=pow(1+1/(double)(k), k) * eta * (x[i]-a) * temp[numXElems*k+i] ; 
	}
	else if(scaling ==0) {
		for(int i=0;i<numXElems*r;i++) temp[i]=1;
		for(int k=0;k<r-1;k++)
			for(int i=0;i<numXElems;i++)
				temp[numXElems*(k+1)+i]= 1/(double)(k+1) * (x[i]-a) * temp[numXElems*k+i] ; 
	}
}

void ComputeB_Scaled(double** outB, int r, double eta0, double a, double b, double dx, double dy, int fun, const int scaling) {

	double ba = b-a;
	//diagS is used in fun values other than 1 and 3. This is a diagonal matrix (rxr) whose anti-diagonal elements is alternating between -1 and 1 starting from 1 at the position diagS(r, 1).
	int tempVal = -1;
	double *diagS=new double[r]; //diagonal matrix stored as a 1D array of r non-zero diagonal elements.
	for(int i=0;i<r;i++) {
		diagS[i]=(i%2)?-1:1;
	}
	
	double *diagDD =new double[r];
	for(int i=0;i<r;i++) {
		diagDD[i]=(i%2)?-1:1;
	}
	
	double etax = eta0 * 2/dx;
	double etay = eta0 * 2/dy;

	*outB = new double[r*r];
	double* B=*outB;
	for(int i=0;i<r*r;i++) B[i]=0;

	if(fun == 1) {
		switch(scaling) {
			case 0:
				{
					for(int k=0;k<r;k++) {
						double temp=-1/(double)ba;
						for(int i=0;i<r-k;i++) {
							B[r*k+i]= temp;
							temp =	temp * (i+1)/ba;
						}
					}
				}
				break;
			case 1:
				{
					double temp=-1/(double)ba;
					B[0]=temp;
					B[1]=temp/(ba*etay);
					B[r]=temp/(ba*etax);
					B[r+1]=2/(ba*etay)*B[r];
					double* temppower = new double[r];
					for(int i=2;i<r;i++)
						temppower[i]=pow(1 - 1/((double)i), i-1);

					//completing first row of B
					for(int i=2;i<r;i++)
						B[i]=1/(double)(ba * etay) * B[i-1] * temppower[i];
					//completing second row of B
					for(int i=2;i<r-1;i++)
						B[r+i]=(i+1)/(double)(ba * etay * i) * B[r+i-1] * temppower[i];

					//completing first column of B
					for(int i=2;i<r;i++)
						B[r*i]=1/(double)(ba * etax) * B[r*(i-1)] * temppower[i];

					//completing second column of B
					for(int i=2;i<r-1;i++)
						B[r*i+1]=(i+1)/(double)(ba * etax * i) * B[r*(i-1)+1] * temppower[i];

					//completing rest of B
					for(int k = 2;k<r-2;k++)
						for(int i = 2;i<r-k+1;i++)
							B[k*r+i] = (k+i)/(double)(ba*etay *i)*B[k*r+i-1]*temppower[i];
					delete [] temppower;
				}
				break;
			default:assert(0);
		}
		
		//Doing B=B*DD;
		for(int i=0;i<r;i++) {
			for(int j=0;j<r;j++) {
				if(j%2)
					B[r*i+j] *= -1;
			}
		}
			
	}
	else if(fun == 2) {
		switch(scaling){
			case 0:
			    {
				double temp=1/(double)ba;
				//B is an rxr array with each diagonal starting from main diagonal has constant k!*-1/ba^k except for k=0
				//creates a symmetric toeplitz square matrix that is triangular above main anti-diagonal. Needed for doing triu(toeplitz(cc(r,r-1,r-2,...,1))) and then multiplying the upper triangular matrix with diagonal matrix whose anti-diagonal elements are non-zeros. B=C*S, where S is an rxr matrix whose only antidiagonal elements are 1 or -1 results in a B, where the columns of C appear in reverse order and multiplied with -1 or 1.
				int fact=1;
				double newDiag=temp;

				//The top left corner is 1/ba^2. The anti-diagonals below the top-left corner element have values +-2!/ba^3, +-3!/ba^4, +-4!/ba^5 and so on.
				for(int i=0;i<r;i++) {
					newDiag=newDiag*temp*fact; 
					for(int j=i;j>=0;j--) {
						B[r*j+(i-j)] = ((i-j)%2)?-newDiag:newDiag;
					}
					fact++;
				}
			    }
				break;
			case 1:
			    {
				double temp=1/(double)ba;
				temp = temp / ba; //i.e. temp is 1/ba^2
				B[0]=temp;
				B[1]=(2/(double)(ba*etay))*temp;
				B[r]=(2/(double)(ba*etax))*temp;
				B[r+1]=3/(ba*etay)*B[r];
				double* temppower = new double[r];
				for(int i=2;i<r;i++)
					temppower[i]=pow(1 - 1/((double)i), i-1);

				//completing first row of B
				for(int i=2;i<r;i++)
					B[i]=(i+1)/(double)(ba * etay * i) * B[i-1] * temppower[i];
				//completing second row of B
				for(int i=2;i<r-1;i++)
					B[r+i]=(i+2)/(double)(ba * etay * i) * B[r+i-1] * temppower[i];

				//completing first column of B
				for(int i=2;i<r;i++)
					B[r*i]=(i+1)/(double)(ba * etax * i)* B[r*(i-1)] * temppower[i];

				//completing second column of B
				for(int i=2;i<r-1;i++)
					B[r*i+1]=(i+2)/(double)(ba * etax * i) * B[r*(i-1)+1] * temppower[i];

				//completing rest of B
				for(int k = 2;k<r-2;k++)
					for(int i = 2;i<r-k+1;i++)
						B[k*r+i] = (k+i)/(double)(ba*etay*i)*B[k*r+i-1]*temppower[i];
				//Doing B=B*DD;
				for(int i=0;i<r;i++) {
					for(int j=0;j<r;j++) {
						if(j%2)
							B[r*i+j] *= -1;
					}
				}
				
				delete [] temppower;

			 }

			 break;

		default: assert(0);
		}
	}
	else if(fun == 3){
		switch(scaling) {
			case 0:
				{
				double temp=-1/(double)ba;
				//cc[r] is an rx1 array containing k!*-1/ba^k
				//creates a symmetric matrix that is triangular above anti diagonal. Needed for doing triu(toeplitz(cc(r,r-1,r-2,...,1)))
				//First create the toeplitxmatrix cc(r, r-1, r-2, ..1), where cc[k]=(-1) * (temp)^(k-1) * (k-2)! where 1<k 
				for(int k=0;k<r;k++) {
					B[k]=temp;
					for(int i=0,j=1;i<k;i++,j++)
						B[r*j+(k-(i+1))]=temp;
					temp=temp * (k+1)/ba;
				}
					
				//adjusting for multiplying with S
				for(int k=0;k<r;k++) 
					for(int i=0;i<r-k;i++)
						if(i % 2) 	
							B[r*k+i] *= -1;

				B[0]=log(abs(b-a));
				}
				break;
			case 1:
				{
				double temp=log(abs(a-b));
				B[0]=temp;
				B[1]=-1/(double)(ba*etay);
				B[r]=-1/(double)(ba*etax);
				B[r+1]=1/(ba*etay)*B[r];
				double* temppower = new double[r];
				for(int i=2;i<r;i++)
					temppower[i]=pow(1 - 1/((double)i), i-1);

				//completing first row of B
				for(int i=2;i<r;i++)
					B[i]=(i-1)/(double)(ba * etay * i) * B[i-1] * temppower[i];
				//completing second row of B
				for(int i=2;i<r-1;i++)
					B[r+i]=1/(double)(ba * etay) * B[r+i-1] * temppower[i];

				//completing first column of B
				for(int i=2;i<r;i++)
					B[r*i]=(i-1)/(double)(ba * etax * i)* B[r*(i-1)] * temppower[i];

				//completing second column of B
				for(int i=2;i<r-1;i++)
					B[r*i+1]=1/(double)(ba * etax) * B[r*(i-1)+1] * temppower[i];

				//completing rest of B
				for(int k = 2;k<r-2;k++)
					for(int i = 2;i<r-k+1;i++)
						B[k*r+i] = (k+i-2)/(double)(ba*etay*i)*B[k*r+i-1]*temppower[i];
				//Doing B=B*DD;
				for(int i=0;i<r;i++) {
					for(int j=0;j<r;j++) {
						if(j%2)
							B[r*i+j] *= -1;
					}
				}
				
				delete [] temppower;
				}
				break;
			default:assert(0);
		}
	}
	else
		printf("TODO: FUN %d NOT SUPPORTED\n",fun);

}


//to compute a function (outer shift) considering two cells. the centers of the cells are a and b. diameters are dx and dy.
void ComputeT_Scaled(double** outT, int r, double eta0, double a, double b, double dx, double dy, const int scaling) {

	*outT = new double[r*r];
	double* T=*outT;

	double temp = dx/dy;
	for(int i=0;i<r*r;i++)
		T[i]=0;

	for(int i=0;i<r;i++)
		T[r*i+i]=pow(temp, i);

	double ab=a-b;
	T[1] = ab*eta0*2/dy;
	
	//assert(abs(T[1]-1*eta0*(dx/2)*ab) == 0); 
		
	switch(scaling) {
		case 0:
			{
			temp=1;
			//creates an upper triangular matrix that is toeplitz with constant ab^k/k! for k=0 to r-1
			for(int k=0;k<r;k++) {
				T[k]=temp;
				for(int i=0,j=1;i<r-k;i++,j++)
					T[r*j+j]=temp;
				temp=temp * ab/k;
			}
			}
			break;

		case 1:
			{
			double* temppower = new double[r];
			for(int i=2;i<r;i++)
				temppower[i]=pow(1 - 1/((double)i), i-1);

			//completing first row of B
			for(int i=2;i<r;i++)
				T[i]=ab * (eta0*2/dy) * T[i-1] / temppower[i];

			//completing rest of B
			for(int k = 1;k<r;k++)
				for(int i = k+1;i<r;i++)
					T[k*r+i] = (ab)/(i-k)*i*(eta0*2/dy)*T[k*r+i-1]/temppower[i];

			delete [] temppower;
			}
			break;
		case 2:
			{
			temp=1;
			//creates an upper triangular matrix that is toeplitz with constant ab^k/k! for k=0 to r-1
			for(int k=0;k<r;k++) {
				T[k]=temp;
				for(int i=0,j=1;i<r-k;i++,j++)
					T[r*j+j]=temp;
				temp=temp * ab/k;
			}
			
			//premultiplying T by diag(1, etax, etax^2, etax^3,...,etax^r-1). This is same as scaling the rows of T by the vector element at the corresponding row index.
			temp=1;
			double scalingFactor=eta0*2/dx;
			for(int i=0;i<r;i++) {
				for(int j=i;j<r;j++)
					T[r*i+j] *= temp;
				temp *= scalingFactor; 
			}

			//postmultiplying T by diag(1, etay, etay^2, etay^3,...,etay^r-1) This is same as scaling the columns of T by the vector element at the corresponding column index.
			temp=1;
			scalingFactor=eta0*2/dy;
			for(int i=0;i<r;i++) {
				for(int j=0;j<=i;j++) 
					T[r*j+i] *= temp;
				temp *= scalingFactor;
			}
			}
			break;
		default:assert(0);
	}
}
