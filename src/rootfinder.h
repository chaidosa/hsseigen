#ifndef ROOTFINDER_H
#define ROOTFINDER_H
#include<vector>
class Root{
    public:
        std::vector<double> tau;
        double *org;
        int org_size;
        std::vector<double> x;
        double percent;
};
Root *rootfinder(std::vector<double>& d, std::vector<double>& v);
#endif
