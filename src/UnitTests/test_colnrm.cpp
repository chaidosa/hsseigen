#include<bits/stdc++.h>
#include<fstream>
#include<iomanip>
#include<algorithm>
#include <iostream>
#include<string.h>
#include "../colnorms.h"
#include "../bsxfun.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}

int main(){
    #ifndef INPUT_OUT
	    freopen("input_colnrm.txt","r",stdin);	
    #endif
    int  n;
    std::cin >> n;
     std::vector<double>d,lam, tau, org, v;
    
    for(int i = 0; i < n; i++){
        double temp;
        std::cin >> temp;
        d.push_back(temp);
    }    

    for(int i = 0; i < n; i++){
        double temp;
        std::cin >> temp;
        lam.push_back(temp);
    }

    for(int i = 0; i < n; i++){
        double temp;
        std::cin >> temp;
        tau.push_back(temp);
    }

    for(int i = 0; i < n; i++){
        double temp;
        std::cin >> temp;
        org.push_back(temp-1);
    }

    for(int i = 0; i < n; i++){
        double temp;
        std::cin >> temp;
        v.push_back(temp);
    }

    std::vector<double>s = colnorms(d, lam, tau, org, v, 1024);
    
    std::ofstream txtOut;
    txtOut.open("output.txt", std::ofstream::out);
    for(int i = 0; i < s.size(); i++)
        txtOut<<s[i]<<std::endl;   

return 0;

}