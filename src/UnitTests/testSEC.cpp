#include<bits/stdc++.h>
#include<fstream>
#include<iomanip>
#include<algorithm>
#include <iostream>
#include<fstream>
#include<iomanip>
#include<string.h>
#include "../secular.h"
#include "../eigenmatrix.h"
using namespace std;
int main(){
    #ifndef INPUT_OUT
	    freopen("/home/pritesh/Music/TEst/SuperDC-main/superdc_1.0.0/tests/input_sec.txt","r",stdin);	
    #endif
    int n = 64;
    
    double *d = new double[n];
    double *v = new double[n];
    for(int i = 0; i < n; i++)
        cin>>d[i];
    
    for(int i = 0; i < n; i++)
        cin>>v[i];
    
    SECU *ret = secular(d, n, v, n, 1024);

    cout<<"Test done successfully";

    return 0;
}