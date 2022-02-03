#include<bits/stdc++.h>
#include<fstream>
#include<iomanip>
#include<algorithm>
#include <iostream>
#include<fstream>
#include<iomanip>
#include<string.h>
#include "../rootfinder.h"
int main(){
    #ifndef INPUT_OUT
	    freopen("input_rootf.txt","r",stdin);	
    #endif
    int  n;
    std::cin >> n;
    std::vector<double>d,v;
    
    for(unsigned int i = 0; i < n; i++){
        double temp;
        std::cin >> temp;
        d.push_back(temp);
    }    

    for(unsigned int i = 0; i < n; i++){
        double temp;
        std::cin >> temp;
        v.push_back(temp);
    }

    Root *r = rootfinder(d,v);
    std::ofstream txtOut;
    txtOut.open("output.txt", std::ofstream::out);
    txtOut<<"x:\n";
    for(unsigned int i = 0; i < r->x.size(); i++)
        txtOut<<r->x[i]<<std::endl;

    txtOut<<"tau:\n";
    for(unsigned int i = 0; i < r->tau.size(); i++)
        txtOut<<r->tau[i]<<std::endl;

    txtOut<<"org:\n";
    for(unsigned int i = 0; i < r->org.size(); i++)
        txtOut<<r->org[i]<<std::endl; 

    txtOut<<"Percent :"<<r->percent;

   
    return 0;
}
