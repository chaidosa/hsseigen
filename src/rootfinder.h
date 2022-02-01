#ifndef ROOTFINDER_H
#define ROOTFINDER_H
#include<vector>
class Root{
    public:
        std::vector<double> tau;
        std::vector<double> org;
        std::vector<double> x;
        double percent;
};
Root *rootfinder(std::vector<double>& d, std::vector<double>& v);
#endif