#include<bits/stdc++.h>
#include "../makeband.h"
using namespace std;
char* testFile;


int main(int argc, char* argv[]){
    char inputFile[]="testinput.txt";
    testFile=inputFile;    
    double* bandMatrix;
    
    //calling this function creates a file testinput.txt and populates it with the band matrix. Also returns the band matrix in the variable bandMatrix
    MakeBand(128, 1, &bandMatrix);
    int n=128;
    for (int i = 0; i < n; i++) {
	    for (int j = 0; j < n; j++) {
			cout<<bandMatrix[i*128+j];
			if(j!=n-1)
				cout<<"\t";
	    }
	    if(i!=n-1)
		cout<<"\n";
  	}
    return 0;
}
