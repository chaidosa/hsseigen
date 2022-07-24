#include<stdio.h>
#include<fstream>
#include<iostream>
#include<math.h>
const double pi = 3.14159265358979323846;
double ba=2.2074e+03; //set the test value here.
int r=50;
double dx=1.7659e+03;
double dy=882.9454;

void test_fun1_case1();
void test_fun2_case0();
void test_fun2_case1();
void test_fun3_case1();

// This for testing the output of matrix
void printArray(double **Arr, int row, int col, const char* filename="output.txt")
{
	std::ofstream txtOut;
    txtOut.open(filename, std::ofstream::out | std::ofstream::app);
    double *A = *Arr;
    for(unsigned int i = 0; i <(unsigned)row; i++)
    {
        for(int j=0; j < col; j++)
        {
            txtOut<<A[j+i*col]<<"\n";
        }
    }
    txtOut.close();
}


int main(int argc, char*argv[]){
		
	//TODO: write code to match expected values and indicate if the test passed. Currently, this is done manually by running through gdb
	//test_fun2_case0();
	//test_fun2_case1();
	//test_fun1_case1();
	test_fun3_case1();
	
	return 0;

}

void test_fun3_case1(){

	double* B=new double[r*r];
	double eta0 = pow(2*pi*r, 0.5/r) / exp(1);
	double etax = eta0 * 2/dx;
	double etay = eta0 * 2/dy;
	
	double temp=log(ba);

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
	
	//printArray(&B, 1, 2500, "computeBoutput.txt");
	delete [] temppower;
}

void test_fun1_case1(){
	double* B=new double[r*r];
	double eta0 = pow(2*pi*r, 0.5/r) / exp(1);
	double etax = eta0 * 2/dx;
	double etay = eta0 * 2/dy;

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


	//Doing B=B*DD;
	for(int i=0;i<r;i++) {
		for(int j=0;j<r;j++) {
			if(j%2)
				B[r*i+j] *= -1;
		}
	}
	delete [] temppower;

}

void test_fun2_case0(){
	
	double temp=1/(double)ba;
	int fact=1;
	double newDiag=temp;
	double* B=new double[r*r];

	for(int i=0;i<r;i++) {
		newDiag=newDiag*temp*fact; 
		for(int j=i;j>=0;j--) {
			//printf("Setting B(%d %d)=B[%d]=%f\n",j,(i-j),r*j+i-j,newDiag);
			B[r*j+(i-j)] = ((i-j)%2)?-newDiag:newDiag;
		}
		fact++;
	}
}

void test_fun2_case1(){

	double temp=1/(double)ba;
	temp = temp / ba; //i.e. temp is 1/ba^2
	double* B=new double[r*r];
	double eta0 = pow(2*pi*r, 0.5/r) / exp(1);
	double etax = eta0 * 2/dx;
	double etay = eta0 * 2/dy;

	
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
		for(int i=2;i<r-k+1;i++)
			B[k*r+i] = (k+i)/(double)(ba*etay*i)*B[k*r+i-1]*temppower[i];

	//Doing B=B*DD;
	for(int i=0;i<r;i++) {
		for(int j=0;j<r;j++) {
			if(j%2)
				B[r*i+j] *= -1;
		}
	}

	//printArray(&B, 1, 2500, "computeBoutput.txt");
	
	delete [] temppower;

}
