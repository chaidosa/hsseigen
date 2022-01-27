#include<bits/stdc++.h>
#include<fstream>
#include<iomanip>
#include<algorithm>
#include <iostream>
#include<fstream>
#include<iomanip>
#include<string.h>
#include "../vhat.h"

int main(){
    #ifndef INPUT_OUT
	    freopen("input_vhat.txt","r",stdin);	
    #endif
    int  n;
    std::cin >> n;
    std::vector<double>d,lam, tau, org, w;
    
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
        w.push_back(temp);
    }

    std::vector<double>v = vhat(d, lam, org, tau, w, 1024);
    
    std::ofstream txtOut;
    txtOut.open("output.txt", std::ofstream::out);
    for(int i = 0; i < v.size(); i++)
        txtOut<<v[i]<<std::endl;   

return 0;
}