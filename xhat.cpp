#include<iostream>
#include<fstream>
#include<cmath>

int main(int argc, char **argv){
    
    std::ifstream xMatlab("X_node_correct.txt");
    std::ifstream xCpp("X_Node.txt");

    bool isPrecise = true;
    double greatest_err = 0.; 

    double rs = 0.;
    int n = 0;
    while ( !xMatlab.eof() && !xCpp.eof() )
    {
        double a = 0;
        double b = 0;

        xCpp >> a;
        xMatlab >> b;

        double error;

        
        error = (a-b)*(a-b);

        //std::cout << "Error is " << error << "\n";

        if ( error > 10e-9)
        {
            isPrecise = false;
        }

        n+=1;
        rs += error;

        error = sqrt(error);
        std::cout << "xCpp is: " << a << " xMatlab is " << b << "\n";
    }
    
    double rmserror = sqrt(rs)/n;

    if ( rmserror > 10e-10){
            isPrecise = false;
    }

    std::cout << "Is the cpp code precise? : " << isPrecise << "\n";
    std::cout << "Rms error is " << rmserror << std::endl;

    xCpp.close();
    xMatlab.close();
    return 0;
}