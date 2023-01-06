#include<iostream>
#include<fstream>
#include<cmath>

int main(int argc, char **argv){
    
    std::ifstream dMatlab("d.txt");
    std::ifstream dCpp("output.txt");

    bool isPrecise = true;

    double rs = 0.;
    int n = 0;
    while ( !dMatlab.eof() && !dCpp.eof() )
    {
        double a = 0;
        double b = 0;

        dCpp >> a;
        dMatlab >> b;

        double error;

        
        error = (a-b)*(a-b);

        //std::cout << "Error is " << error << "\n";

        if ( error > 10e-15)
        {
            isPrecise = false;
        }

        n+=1;
        rs += error;
    }
    
    double rmserror = sqrt(rs)/n;

    std::cout << "Is the cpp code precise? : " << isPrecise << "\n";
    std::cout << "Rms error is " << rmserror << std::endl;

    dCpp.close();
    dMatlab.close();
    return 0;
}
