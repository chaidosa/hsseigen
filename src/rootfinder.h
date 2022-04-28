#ifndef ROOTFINDER_H
#define ROOTFINDER_H
#include<vector>
class Root{
    public:
        std::vector<double> tau;
        int *org;
        int org_size;
        double* x;
        double percent;
	Root(int p_OrgSize) {
		org_size=p_OrgSize;
		x=NULL;
		org=NULL;
	}
};
Root *rootfinder(std::vector<double>& d, std::vector<double>& v);
#endif
