#include<iostream>
#include<string.h>
#include<stdlib.h>
#include<fstream>
#include"RandGen.h"
#include"makeband.h"
//#include<random>

using namespace std;
extern char* testFile;
#ifdef DEBUGINPUT
int counter=1.; 
#endif
/* Function to create a symmetric banded matrix. e.g. creating a matrix of bandwidth 3 (r) and 6 (n) rows
 * 1 2 3  0  0  0
 * 2 4 5  6  0  0
 * 3 5 7  8  9  0
 * 0 6 8 10 11 12
 * 0 0 9 11 13 14
 * 0 0 0 12 14 15
 *
 * the created matrix is allocated memory and returned in bandMatrix.
 * if RANDOM_GEN is defined the values are initialized with random values. Otherwise, the values are read
 * from a file called test_input.txt. 
 * Note that for debug purposes, there is a special initialization sequence.
 */ 
int MakeBand(int n, int r, double** bandMatrix)
{
  int i;
  int j;
  //double b_j;
  *bandMatrix = new double[n*n];
  double* A = *bandMatrix;
  memset(A, 0, sizeof(double)*n*n);

  //using C++11 features	
  /*std::random_device rd;  //Will be used to obtain a seed for the random number engine
  //std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 gen;
  gen.seed(0);
  std::uniform_real_distribution<> dis(0.0, 1.0);*/
#ifdef RANDOM_GEN
  //The default generator of MATLAB is mersenne twister with seed 0.
  //Using wrapper for mt19937 C++11 not required.
  RandGen myrand(UNIFORM_REAL);
  myrand.SetSeed(0);
  myrand.SetInterval(0.0,1.0);

  for (i = 0; i < n; i++) {
    int endJ = std::min(i+r,n);
    for (j = i; j < endJ; j++) {
        	A[j + i * n] = 2 * (myrand.Rand() - 0.5);
        	//A[j + i * n] = 2 * (dis(gen) - 0.5);
        	A[i + j * n] = A[j + i*n];
    }
  }
#elif DEBUGINPUT
  std::ofstream outputStr(testFile,std::ofstream::out);
  if(!outputStr.is_open()) {
	cout<<"Error opening"<<testFile<<endl;
	return -1;
  }


  for (i = 0; i < n; i++) {
    for(int ii=0;ii<i;ii++)
	    outputStr<<"0\t";
    int endJ = std::min(i+r,n);
    for (j = i; j < endJ; j++) {
        	A[j + i * n] = counter++;
        	A[i + j * n] = A[j + i*n];
		//outputStr<<"["<<i<<","<<j<<"]"<<A[i+j*n];
		outputStr<<A[i+j*n];
		if(j!=(n-1))
			outputStr<<"\t";
    }
    for(int ii=endJ;ii<n;ii++)
	    outputStr<<"0\t";
    if(i!=(n-1))
    	outputStr<<"\n";
  }
  outputStr.close();
#else
  std::ifstream inputStr(testFile,std::ifstream::in);

  if(!inputStr.is_open()) {
	cout<<"Error opening test_input.txt."<<endl;
	return -1;
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
		double val;
		inputStr>>val;
        	A[i*n + j] = val;
       // 	A[i + j * n] = A[j + i*n];
    }
  }
  inputStr.close();
#endif
  return 0;
}

