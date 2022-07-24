#include<cstdlib>
#include<cassert>
#include<cmath>
#include<set>
#include<fstream>
#include<iomanip>
#include"computeX_Scaled.h"
#include"fmm_types.h"
#include "QR.h"
#include "bsxfun.h"
#ifndef OPENBLAS
extern "C"
{
#endif
#include<cblas.h>
#ifndef OPENBLAS
}
#endif

const double TAU=0.6; //also called separation ratio or opening radius (in case of tree representing 2D/3D points) 
const double pi = 3.14159265358979323846;

using namespace std;

void PrintTree(const Vertex* node, std::ofstream& of);
void DestroyTree(Vertex* node);
class TreeVisitor{

	public:
	Vertex* rootNode;
	double eta0;
	const double* dataX, *dataY; //These are local copies of arguments (X and Y vectors) passed to the fmm_xx functions. Used within this module/file only.
	int scalingLocal, numTerms, funLocal; //These are local copies of arguments passed to the fmm_xx functions. Used within this module/file only.
	const double* dataQ;//This is used within this module/file to refer to the memory locations that store the Q (e.g. charge) array
	const int* org; //This is an array of indices that indicate the shift to be used.
	const double* gap; //This is used within this module/file to refer to the memory locations that indicate difference between successive elemnts (X[i]-Y[i]) or (Y[i]-X[i]) or .. 
	const int orgSize; //Size of the org array. also indicates the number of elements in the X vector. This is used to determine the size of the result. 

	//'result' points to the result computed by fmm functions. Always of size numXElems (a parameter to the fmm function). 
	//In case of triangular fmm function, upper half is appended to lower half. So, the size is 2* numXElems.
	double* result; 
	
	TreeVisitor(const double* x, const double* y, const double* q, int r, const int scaling, const int fun, const int * p_org, const double* p_gap, const int p_orgSize, Vertex* node):dataX(x), dataY(y), dataQ(q), numTerms(r), scalingLocal(scaling), funLocal(fun), org(p_org), gap(p_gap), orgSize(p_orgSize), rootNode(node){
		eta0 = pow(2*pi*r, 0.5/r) / exp(1);
		result=NULL;
	}
	virtual void LeafNodeActions(Vertex* node);
	virtual void PreorderActions(Vertex* node);
	virtual void PostorderActions(Vertex* node);
	virtual ~TreeVisitor();
};

/*
template<typename T>
void PrintArray(T *A, int row, int col, const char* filename="debugoutput.txt")
{
#if 1
    ofstream txtOut;
    txtOut.open(filename, std::ofstream::out | std::ofstream::app);
    for(unsigned int i = 0; i <(unsigned)row; i++)
    {
        for(int j=0; j < col; j++)
        {
            txtOut<<std::setprecision(12)<<A[j+i*col]<<"\n";
        }
    }
    txtOut.close();
#endif
}

template<typename T>
void PrintArray2(T *A, int row, int col, const char* filename="debugoutput.txt")
{
#if 1
    ofstream txtOut;
    txtOut.open(filename, std::ofstream::out | std::ofstream::app);
    for(unsigned int i = 0; i <(unsigned)row; i++)
    {
        for(int j=0; j < col; j++)
        {
            txtOut<<std::setprecision(12)<<A[j+i*col];
		if(j != col-1)
			txtOut<<" ";
        }
	if(i != row-1)
		txtOut<<"\n";
    }
    txtOut.close();
#endif
}

*/

void TreeVisitor::LeafNodeActions(Vertex* node) {
	if((node->left==NULL) && (node->right==NULL)) {
		if((node->xRight - node->xLeft) != 0) {
			double* outU=NULL;
			int numElements=node->xRight-node->xLeft;
			if(numElements > 0) {
				ComputeU_Scaled(&outU, dataX+node->xLeft, (node->xRight-node->xLeft), numTerms, eta0, node->center, 2*(node->radius), scalingLocal);
				GetTransposeInPlace(outU, numTerms, (node->xRight-node->xLeft));
				if(node->computed) {
					//compute z(px, :) = Ui * u{i};
					node->zl=new double[numElements];
					for(int j=0;j<numElements;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += outU[j*numTerms+i] * node->ul[i]; 
						}
						node->zl[j] = temp; 
						assert((node->xLeft+j) < orgSize);
						result[node->xLeft+j]=temp;
					}
				}
				delete [] outU;
			} else{
				assert(0);
			}
			//Union(i, neighbor{i})
			bool isSelf=true; 
			node->neighbors.insert(node->neighbors.begin(), node);
			for(std::vector<Vertex*>::iterator iter=node->neighbors.begin();iter!=node->neighbors.end();iter++) {
				Vertex* nbr = *iter;
				int xjVectorSize = nbr->yRight-nbr->yLeft;
				if(xjVectorSize !=0) {
					double* tj = new double[xjVectorSize];
					double* xj = new double[xjVectorSize];
					for(int i=0;i<xjVectorSize;i++) {
						assert(i < orgSize);
						xj[i] = dataX[org[nbr->yLeft+i]]; 
						tj[i] = gap[nbr->yLeft+i];
					}
					int xiVectorSize = node->xRight-node->xLeft;
					double *D = new double[xiVectorSize*xjVectorSize];
					//computing D=xi-xj.' 
					for(int i=0;i<xiVectorSize;i++) {
						for(int j=0;j<xjVectorSize;j++) {
					 		D[i*xjVectorSize+j] = dataX[node->xLeft+i] - xj[j];
						}
					}
					delete [] xj;
					bsxfun('m',&D, std::make_pair(xiVectorSize, xjVectorSize), tj, std::make_pair(1,xjVectorSize)); 
					delete [] tj;
					switch(funLocal) {
						case 1: 
							{
								for(int i=0;i<xiVectorSize;i++) 
									for(int j=0;j<xjVectorSize;j++) 
					 					D[i*xjVectorSize+j] = 1/(double)D[i*xjVectorSize+j];
							}
							break;
						case 2:
							{
								for(int i=0;i<xiVectorSize;i++) 
									for(int j=0;j<xjVectorSize;j++) 
					 					D[i*xjVectorSize+j] = 1/(double)(D[i*xjVectorSize+j] * D[i*xjVectorSize+j]);
							}
							break;
						case 3:
							{
								for(int i=0;i<xiVectorSize;i++) { 
									for(int j=0;j<xjVectorSize;j++) {
					 					D[i*xjVectorSize+j] = std::log(std::abs(D[i*xjVectorSize+j]));
										if(isSelf && isinf(D[i*xjVectorSize+j]))
											D[i*xjVectorSize+j] = 0;
									}
								}
							}
							break;
						default:
							assert(0);				
					}

					//cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, xiVectorSize, 1, xjVectorSize, 1, D, xjVectorSize, dataQ+nbr->yLeft, 1, 1, node->res, 1);	
					cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, xiVectorSize, 1, xjVectorSize, 1, D, xjVectorSize, dataQ+nbr->yLeft, 1, 1, result+node->xLeft, 1);	
                    			//z(px, :) = z(px, :) + D * q(py, :);
					delete [] D;
				}
				isSelf=false; //set flag to false in all iterations except the first one. When this is false, (isnan() call is not done)
			}
			node->neighbors.erase(node->neighbors.begin());
		}
	}
	
}

void TreeVisitor::PreorderActions(Vertex* node) {
	if(node->level >= 2) {
		for(std::vector<Vertex*>::iterator iter=node->wellSeparatedNodes.begin();iter!=node->wellSeparatedNodes.end();iter++) {
			double* outB; 
			Vertex* wsNode=*iter;
			//printf("ComputeB called with eta0:%lf a:%lf b:%lf, dx:%lf, dy:%lf, fun:%d, scaling:%d\n",eta0, node->center, wsNode->center, 2*(node->radius),2*(wsNode->radius),funLocal, scalingLocal); 
			/*PrintArray<double>(&eta0, 1, 1, "fmminput_Bscaled.txt");
			PrintArray<double>(&(node->center), 1, 1, "fmminput_Bscaled.txt");
			PrintArray<double>(&(wsNode->center), 1, 1, "fmminput_Bscaled.txt");
			double arg1=2*(node->radius);
			double arg2=2*(wsNode->radius);
			PrintArray<double>(&arg1, 1, 1, "fmminput_Bscaled.txt");
			PrintArray<double>(&arg2, 1, 1, "fmminput_Bscaled.txt");
			PrintArray<int>(&funLocal, 1, 1, "fmminput_Bscaled.txt");
			PrintArray<int>(&scalingLocal, 1, 1, "fmminput_Bscaled.txt");*/
            		ComputeB_Scaled(&outB, numTerms, eta0, node->center, wsNode->center, 2*(node->radius), 2*(wsNode->radius), funLocal, scalingLocal);
			//PrintArray2<double>(outB, numTerms, numTerms, "fmmoutput_Bscaled.txt");
			if(node->computed == false) {
				assert(wsNode->res);
				//computing u{i}=Bij * v{j}
				node->ul=new double[numTerms];
				for(int j=0;j<numTerms;j++){
					double temp=0;
					for(int i=0;i<numTerms;i++){
						temp += (outB[j*numTerms+i] * wsNode->res[i]); 
					}
					node->ul[j] = temp;
				}
				node->computed=true;
			} else {
				//computing u{i}+=Bij * v{j}
				for(int j=0;j<numTerms;j++){
					double temp=0;
					for(int i=0;i<numTerms;i++){
						temp += (outB[j*numTerms+i] * wsNode->res[i]); 
					}
					node->ul[j] += temp;
				}
			}
			delete [] outB;
		}
	

		if(node->level > 3) {
		    Vertex* parent = node->parent;
		    if(parent->computed) {
			assert(parent->ul != NULL);
			double *outT=NULL;
			ComputeT_Scaled(&outT, numTerms, eta0, node->center, parent->center,  2*(node->radius), 2*(parent->radius), scalingLocal);
			if(node->computed == false) {
				assert(node->ul == NULL);
				//computing u{i}=Ri * u{p}
				node->ul=new double[numTerms];
				for(int j=0;j<numTerms;j++){
					double temp=0;
					for(int i=0;i<numTerms;i++){
						temp += (outT[j*numTerms+i] * parent->ul[i]); 
					}
					node->ul[j] = temp;
				}
				node->computed=true;
			} else {
				//computing u{i}+=Ri * u{p}
				for(int j=0;j<numTerms;j++){
					double temp=0;
					for(int i=0;i<numTerms;i++){
						temp += (outT[j*numTerms+i] * parent->ul[i]); 
					}
					node->ul[j] += temp;
				}
			}
			delete [] outT;
		    }
		}	    
	}
}

void TreeVisitor::PostorderActions(Vertex* node) {
    if(node->level >= 2) {
        if((node->left==NULL) && (node->right==NULL)) {
		double* yNode=NULL;
		//allocates memory and populates yNode (dimension: (yRb-yLb) x numTerms i.e yNode dimensions are transposed) 
		
		int numElements=node->yRight-node->yLeft;
		if(numElements > 0) {
			ComputeU_Scaled(&yNode, dataY+node->yLeft, numElements, numTerms, eta0, node->center, 2*(node->radius), scalingLocal);
			//computing matrix-vector prod of yNode and q to store the result in yNode)
			node->res = new double[numTerms];
			for(int j=0;j<numTerms;j++){
				double temp=0;
				for(int i=0;i<numElements;i++){
					temp += yNode[j*numElements+i] * dataQ[node->yLeft+i]; 
				}
				node->res[j] = temp;
			}
			delete [] yNode;
		}
		else {
			node->res = new double[numTerms];
			for(int i=0;i<numTerms;i++)
				node->res[i]=0;
			/*yNode = new double[numTerms* numTerms];
			//vec: 
			for(int i=0;i<numTerms*numTerms;i++)
				yNode[i]=0;
			node->res=yNode; //Now first numTerms of res array contain v{i} */
			//assert(0);
		}
		//node->computed=true;
	} else {
		double* nodeLeftOutput=NULL, *nodeRightOutput=NULL;
	    	assert(node->left && node->right);
	    	ComputeT_Scaled(&nodeLeftOutput, numTerms, eta0, node->left->center, node->center, 2*(node->left->radius), 2*(node->radius), scalingLocal);
	    	ComputeT_Scaled(&nodeRightOutput, numTerms, eta0, node->right->center, node->center, 2*(node->right->radius), 2*(node->radius), scalingLocal);
		GetTransposeInPlace(nodeLeftOutput, numTerms, numTerms);
		GetTransposeInPlace(nodeRightOutput, numTerms, numTerms);
		assert(node->res == NULL);
		node->res=new double[numTerms];
		for(int j=0;j<numTerms;j++){
            		double temp=0;
			for(int i=0;i<numTerms;i++){
				assert(node->left && node->left->res);
				assert(node->right && node->right->res);
				temp += (nodeLeftOutput[j*numTerms+i] * node->left->res[i]) + (nodeRightOutput[j*numTerms+i] * node->right->res[i]); 
			}
			node->res[j] = temp;
		}
		//node->computed=true;
		delete [] nodeRightOutput;
		delete [] nodeLeftOutput;
    	}
    }
}

TreeVisitor::~TreeVisitor(){
	DestroyTree(rootNode);
	delete rootNode;
}

class TriFMM1TreeVisitor:public TreeVisitor{
	public:
		TriFMM1TreeVisitor(const double* x, const double* y, const double* q, int r, int scaling, int fun, const int * p_org, const double* p_gap, int p_orgSize, Vertex* node):TreeVisitor(x, y, q, r,scaling, fun, p_org, p_gap, p_orgSize, node){}
		void PreorderActions(Vertex* node);
		void LeafNodeActions(Vertex* node);
		virtual ~TriFMM1TreeVisitor(){}
};

class FMM1TreeVisitor:public TreeVisitor{

	public:
		FMM1TreeVisitor(const double* x, const double* y, const double* q, int r, int scaling, int fun, const int * p_org, const double* p_gap, int p_orgSize, Vertex* node):TreeVisitor(x, y, q, r,scaling, fun, p_org, p_gap, p_orgSize, node){}
		void LeafNodeActions(Vertex* node);
		virtual ~FMM1TreeVisitor(){}
};


int compare(const void* x, const void* y) {
	const double arg1 = *static_cast<const double*>(x);
	const double arg2 = *static_cast<const double*>(y);
	if (arg1 < arg2) return -1;
	if (arg1 > arg2) return 1;
	return 0;
}

void DestroyTree(Vertex* node)
{
	if(node->left){
		DestroyTree(node->left);
		delete node->left;
		node->left=NULL;
	}

	if(node->right) {
		DestroyTree(node->right);
		delete node->right;
		node->right=NULL;
	}

	return;
}

Vertex* ConstructSpatialTree(double *x, double *y, int xLb, int xUb, int yLb, int yUb, double center, double rad, int el, int er, Vertex* parent) {

	double c = center - rad/2; //new center
	double d = rad/2; //new radius
	double mid=0.;
	bool flag=false;
	int i, j, k, l, countXLeft=0, countXRight=0, countYLeft=0, countYRight=0;
	Vertex* leftChild=NULL, *rightChild=NULL, *node=NULL;

	if((xLb > xUb) || (yLb > yUb)) {
		assert(0);
		return NULL;
	}

	int sizex = xUb-xLb, sizey = yUb-yLb; //calculate number of elements in the interval
	if ((sizex <= MAX_POINTS_IN_CELL) && (sizey <= MAX_POINTS_IN_CELL)){
		//if this function is called with arguments such that the number of elements in an interval is below a threshold, then stop subdividing the interval.
		node = new Vertex();
		node->parent = parent;
		node->xLeft = xLb;
		node->xRight = xUb;
		node->yLeft = yLb;
		node->yRight = yUb;
		node->isLeaf = true;
		node->level = parent->level + 1;
		node->center = center; //note: not the new center.
		node->radius = rad; //note: not the new radius
		//node->eta = eta0/rad;
		return node;
	}

	//el is the flag indicating if left boundary is to be included. get count of how many elements end up in the left interval (after division of the current interval) 
	if(el==1){
		for(i=xLb;i<xUb;i++) {
			if(((center-rad)<=x[i]) && (x[i]<=center))
				countXLeft++;
		}
		for(j=yLb;j<yUb;j++) {
			if(((center-rad)<=y[j]) && (y[j]<=center))
				countYLeft++;

		}
	} else if (el == 0) {
		for(i=xLb;i<xUb;i++) {
			if(((center-rad)<x[i]) && (x[i]<=center))
				countXLeft++;
		}
		for(j=yLb;j<yUb;j++){
			if(((center-rad)<y[j]) && (y[j]<=center))
				countYLeft++;
		}
	}

	//el is the flag indicating if left boundary is to be included. get count of how many elements end up in the right interval (after division of the current interval) 
	if(er==1) {
		for(k=xLb;k<xUb;k++)
			if((center<x[k]) && (x[k]<=(center+rad)))
				countXRight++;
		for(l=yLb;l<yUb;l++)
			if((center<y[l]) && (y[l]<=(center+rad)))
				countYRight++;
	} else if (er == 0) {
		for(k=xLb;k<xUb;k++)
			if(((center<x[k]) && (x[k]<(center+rad))))
				countXRight++;
		for(l=yLb;l<yUb;l++);
			if((center<y[l]) && (y[l]<(center+rad)))
				countYRight++;
	}
	
	if((countXLeft == 0) || (countYLeft == 0) || (countXRight == 0) || (countYRight == 0)){
		countXLeft=0; countXRight=0; countYLeft=0; countYRight=0;
		//x and y intervals are empty after the bisected point. Adjust center and rad.
		if(((xUb-xLb) > 0) && ((yUb-yLb) > 0))
			mid = (x[xUb-1]+y[yLb])/2;
		else if ((xUb-xLb) > 0)
			mid = (center+rad+x[xLb]) / 2;
		else if  ((yUb-yLb) > 0)
			mid = (center-rad+y[yUb-1]) / 2;
		c = (mid + (center - rad))/2;
		d = (mid - (center - rad))/2;
		if(el==1){
			for(i=xLb;i<xUb;i++)
				if(((center-rad)<=x[i]) && (x[i]<=mid))
					countXLeft++;
			for(j=yLb;j<yUb;j++)
				if(((center-rad)<=y[j]) && (y[j]<=mid))
					countYLeft++;
		} else if (el == 0) {
			for(i=xLb;i<xUb;i++)
				if(((center-rad)<x[i]) && (x[i]<=mid))
					countXLeft++;
			for(j=yLb;j<yUb;j++)
				if(((center-rad)<y[j]) && (y[j]<=mid))
					countYLeft++;
		}

		if(er==1) {
			for(k=xLb;k<xUb;k++)
				if((mid<x[k]) && (x[k]<=(center+rad)))
					countXRight++;
			for(l=yLb;l<yUb;l++)
				if((mid<y[l]) && (y[l]<=(center+rad)))
					countYRight++;
		} else if (er == 0) {
			for(k=xLb;k<xUb;k++)
				if((mid<x[k]) && (x[k]<(center+rad)))
					countXRight++;
			for(l=yLb;l<yUb;l++)
				if((mid<y[l]) && (y[l]<(center+rad)))
					countYRight++;
		}
		flag=true;
	}

	node = new Vertex();
	node->parent = parent;
	node->level = parent->level + 1;
	node->center = center;
	node->radius = rad;
	
	assert((xUb-countXRight) == (xLb+countXLeft));
	assert((yUb-countYRight) == (yLb+countYLeft));

	if((countXLeft > MAX_POINTS_IN_CELL) || (countYLeft > MAX_POINTS_IN_CELL)) {
		node->left = ConstructSpatialTree(x,y,xLb,xLb+countXLeft,yLb,yLb+countYLeft, c, d, el, 1, node); 
	} else{
		node->left = new Vertex();
		node->left->parent = node;
		node->left->xLeft = xLb;
		node->left->xRight = xLb+countXLeft;
		node->left->yLeft = yLb;
		node->left->yRight = yLb+countYLeft;
		node->left->isLeaf = true;
		node->left->level = node->level + 1;
		node->left->center = center; 
		node->left->radius = rad; 
	}

	
	c = center + rad/2; //if(!flag) then the value of c.
	if(flag) {
		c = (mid + (center + rad))/2; //overwrite the value of c if flag is true.
		d = ((center + rad) - mid)/2; 
	}
	if((countXRight > MAX_POINTS_IN_CELL) || (countYRight > MAX_POINTS_IN_CELL)) {
		node->right = ConstructSpatialTree(x,y,xUb-countXRight,xUb,yUb-countYRight,yUb, c, d, 0, er, node); 
	} else{
		node->right = new Vertex();
		node->right->parent = node;
		node->right->xLeft = xLb+countXLeft;
		node->right->xRight = xUb;
		node->right->yLeft = yLb+countYLeft;
		node->right->yRight = yUb;
		node->right->isLeaf = true;
		node->right->level = node->level + 1;
		node->right->center = center; 
		node->right->radius = rad; 
	}

	return node;
}

Vertex* ConstructSpatialTree_trifmm_local_shift(double *x, double *y, int xLb, int xUb, int yLb, int yUb, double center, double rad, int el, int er, Vertex* parent) {

	double c = center - rad/2; //new center
	double d = rad/2; //new radius
	double mid;
	bool flag=false;
	int i, j, k, l, countXLeft=0, countXRight=0, countYLeft=0, countYRight=0;
	Vertex* leftChild=NULL, *rightChild=NULL, *node=NULL;

	if((xLb > xUb) || (yLb > yUb)) {
		assert(0);
		return NULL;
	}

	int sizex = xUb-xLb, sizey = yUb-yLb; //calculate number of elements in the interval
	if ((sizex <= MAX_POINTS_IN_CELL) && (sizey <= MAX_POINTS_IN_CELL)){
		//if this function is called with arguments such that the number of elements in an interval is below a threshold, then stop subdividing the interval.
		node = new Vertex();
		node->parent = parent;
		node->xLeft = xLb;
		node->xRight = xUb;
		node->yLeft = yLb;
		node->yRight = yUb;
		node->isLeaf = true;
		node->level = parent->level + 1;
		node->center = center; //note: not the new center.
		node->radius = rad; //note: not the new radius
		//node->eta = eta0/rad;
		return node;
	}

	//el is the flag indicating if left boundary is to be included. get count of how many elements end up in the left interval (after division of the current interval) 
	if(el==1){
		for(i=xLb;i<xUb;i++) {
			if(((center-rad)<=x[i]) && (x[i]<=center))
				countXLeft++;
		}
		for(j=yLb;j<yUb;j++) {
			if(((center-rad)<=y[j]) && (y[j]<=center))
				countYLeft++;

		}
	} else if (el == 0) {
		for(i=xLb;i<xUb;i++) {
			if(((center-rad)<x[i]) && (x[i]<=center))
				countXLeft++;
		}
		for(j=yLb;j<yUb;j++){
			if(((center-rad)<y[j]) && (y[j]<=center))
				countYLeft++;
		}
	}

	//el is the flag indicating if left boundary is to be included. get count of how many elements end up in the right interval (after division of the current interval) 
	if(er==1) {
		for(k=xLb;k<xUb;k++)
			if((center<x[k]) && (x[k]<=(center+rad)))
				countXRight++;
		for(l=yLb;l<yUb;l++)
			if((center<y[l]) && (y[l]<=(center+rad)))
				countYRight++;
	} else if (er == 0) {
		for(k=xLb;k<xUb;k++)
			if(((center<x[k]) && (x[k]<(center+rad))))
				countXRight++;
		for(l=yLb;l<yUb;l++);
			if((center<y[l]) && (y[l]<(center+rad)))
				countYRight++;
	}
	
	if((countXLeft == 0) || (countYLeft == 0) || (countXRight == 0) || (countYRight == 0)){
		countXLeft=0; countXRight=0; countYLeft=0; countYRight=0;
		//x and y intervals are empty after the bisected point. Adjust center and rad.
		if(((xUb-xLb) > 0) && ((yUb-yLb) > 0))
			mid = (y[yUb-1]+y[yLb])/2;
		else if ((xUb-xLb) > 0)
			mid = (center+rad+x[xLb]) / 2;
		else if  ((yUb-yLb) > 0)
			mid = (center-rad+y[yUb-1]) / 2;
		c = (mid + (center - rad))/2;
		d = (mid - (center - rad))/2;
		if(el==1){
			for(i=xLb;i<xUb;i++)
				if(((center-rad)<=x[i]) && (x[i]<=mid))
					countXLeft++;
			for(j=yLb;j<yUb;j++)
				if(((center-rad)<=y[j]) && (y[j]<=mid))
					countYLeft++;
		} else if (el == 0) {
			for(i=xLb;i<xUb;i++)
				if(((center-rad)<x[i]) && (x[i]<=mid))
					countXLeft++;
			for(j=yLb;j<yUb;j++)
				if(((center-rad)<y[j]) && (y[j]<=mid))
					countYLeft++;
		}

		if(er==1) {
			for(k=xLb;k<xUb;k++)
				if((mid<x[k]) && (x[k]<=(center+rad)))
					countXRight++;
			for(l=yLb;l<yUb;l++)
				if((mid<y[l]) && (y[l]<=(center+rad)))
					countYRight++;
		} else if (er == 0) {
			for(k=xLb;k<xUb;k++)
				if((mid<x[k]) && (x[k]<(center+rad)))
					countXRight++;
			for(l=yLb;l<yUb;l++)
				if((mid<y[l]) && (y[l]<(center+rad)))
					countYRight++;
		}
		flag=true;
	}

#if 1
	if((countXLeft > 0) || (countYLeft > 0)) {
		node = new Vertex();
		node->parent = parent;
		node->level = parent->level + 1;
		/*node->xLeft = xLb;
		node->xRight = xUb;
		node->yLeft = yLb;
		node->yRight = yUb;*/
		leftChild = ConstructSpatialTree_trifmm_local_shift(x,y,xLb,xLb+countXLeft,yLb,yLb+countYLeft, c, d, el, 1, node); 
	}
	else {
		printf("xLb:%d xUb:%d yLb:%d yUb:%d mid:%lf center:%lf radius:%lf\n",xLb, xUb, yLb, yUb, mid, center, rad);
		assert(0); //should not come here.
	}
#endif
/*	node = new Vertex();
	node->parent = parent;
	node->level = parent->level + 1;
	leftChild = ConstructSpatialTree(x,y,xLb,xLb+countXLeft,yLb,yLb+countYLeft, c, d, el, 1, node); */

	c = center + rad/2; //if(!flag) then the value of c.
	if(flag) {
		c = (mid + (center + rad))/2; //overwrite the value of c if flag is true.
		d = ((center + rad) - mid)/2; 
	}
	
	assert(node != NULL);
	rightChild = ConstructSpatialTree_trifmm_local_shift(x,y,xLb+countXLeft,xUb,yLb+countYLeft,yUb, c, d, 0, er, node); 
	node->left = leftChild;
	node->right = rightChild;
	node->center = center;
	node->radius = rad;
	//node->eta = eta0/rad;
	return node;
}
/* x and y contain elements at indices [xLb,xUb) [yLb,yUb). center is the absolute value of coordinate of the center of the interval. rad is the radius if the interval is bisected. */
Vertex* ConstructSpatialTree_temp(double *x, double *y, int xLb, int xUb, int yLb, int yUb, double center, double rad, int el, int er, Vertex* parent) {

	double c = center - rad/2; //new center
	double d = rad/2; //new radius
	double mid;
	bool flag=false;
	int i, j, k, l, countXLeft=0, countXRight=0, countYLeft=0, countYRight=0;
	Vertex* leftChild=NULL, *rightChild=NULL, *node=NULL;

	if((xLb > xUb) || (yLb > yUb)) {
		assert(0);
		return NULL;
	}
	
	std::cout<<setprecision(12)<<"x("<<x[xLb]<<", "<<x[xUb-1]<<") y("<<y[yLb]<<", "<<y[yUb-1]<<")\n";
	int sizex = xUb-xLb, sizey = yUb-yLb; //calculate number of elements in the interval
	if ((sizex <= MAX_POINTS_IN_CELL) && (sizey <= MAX_POINTS_IN_CELL)){
		//if this function is called with arguments such that the number of elements in an interval is below a threshold, then stop subdividing the interval.
		node = new Vertex();
		node->parent = parent;
		node->xLeft = xLb;
		node->xRight = xUb;
		node->yLeft = yLb;
		node->yRight = yUb;
		node->isLeaf = true;
		node->level = parent->level + 1;
		node->center = center; //note: not the new center.
		node->radius = rad; //note: not the new radius
		//node->eta = eta0/rad;
		return node;
	}

	//el is the flag indicating if left boundary is to be included. get count of how many elements end up in the left interval (after division of the current interval) 
	if(el==1){
		for(i=xLb;i<xUb;i++) {
			if(((center-rad)<=x[i]) && (x[i]<=center))
				countXLeft++;
		}
		for(j=yLb;j<yUb;j++) {
			if(((center-rad)<=y[j]) && (y[j]<=center))
				countYLeft++;

		}
	} else if (el == 0) {
		for(i=xLb;i<xUb;i++) {
			if(((center-rad)<x[i]) && (x[i]<=center))
				countXLeft++;
		}
		for(j=yLb;j<yUb;j++){
			if(((center-rad)<y[j]) && (y[j]<=center))
				countYLeft++;
		}
	}

	//el is the flag indicating if left boundary is to be included. get count of how many elements end up in the right interval (after division of the current interval) 
	if(er==1) {
		for(k=xLb;k<xUb;k++)
			if((center<x[k]) && (x[k]<=(center+rad)))
				countXRight++;
		for(l=yLb;l<yUb;l++)
			if((center<y[l]) && (y[l]<=(center+rad)))
				countYRight++;
	} else if (er == 0) {
		for(k=xLb;k<xUb;k++)
			if(((center<x[k]) && (x[k]<(center+rad))))
				countXRight++;
		for(l=yLb;l<yUb;l++);
			if((center<y[l]) && (y[l]<(center+rad)))
				countYRight++;
	}
	
	if((countXLeft == 0) || (countYLeft == 0) || (countXRight == 0) || (countYRight == 0)){
		countXLeft=0; countXRight=0; countYLeft=0; countYRight=0;
		//x and y intervals are empty after the bisected point. Adjust center and rad.
		if(((xUb-xLb) > 0) && ((yUb-yLb) > 0)) {
			mid = (y[yUb-1]+x[xLb])/2;
		}
		else if ((xUb-xLb) > 0){
			mid = (center+rad+x[xLb]) / 2;
			
		}
		else if  ((yUb-yLb) > 0)
			mid = (center-rad+y[yUb-1]) / 2;
		else
			assert(0);
		c = (mid + (center - rad))/2;
		d = (mid - (center - rad))/2;
		if(el==1){
			for(i=xLb;i<xUb;i++)
				if(((center-rad)<=x[i]) && (x[i]<=mid))
					countXLeft++;
			for(j=yLb;j<yUb;j++)
				if(((center-rad)<=y[j]) && (y[j]<=mid))
					countYLeft++;
		} else if (el == 0) {
			for(i=xLb;i<xUb;i++)
				if(((center-rad)<x[i]) && (x[i]<=mid))
					countXLeft++;
			for(j=yLb;j<yUb;j++)
				if(((center-rad)<y[j]) && (y[j]<=mid))
					countYLeft++;
		}

		if(er==1) {
			for(k=xLb;k<xUb;k++)
				if((mid<x[k]) && (x[k]<=(center+rad)))
					countXRight++;
			for(l=yLb;l<yUb;l++)
				if((mid<y[l]) && (y[l]<=(center+rad)))
					countYRight++;
		} else if (er == 0) {
			for(k=xLb;k<xUb;k++)
				if((mid<x[k]) && (x[k]<(center+rad)))
					countXRight++;
			for(l=yLb;l<yUb;l++)
				if((mid<y[l]) && (y[l]<(center+rad)))
					countYRight++;
		}
		flag=true;
	}

#if 1
	if((countXLeft > 0) || (countYLeft > 0)) {
		node = new Vertex();
		node->parent = parent;
		node->level = parent->level + 1;
		/*node->xLeft = xLb;
		node->xRight = xUb;
		node->yLeft = yLb;
		node->yRight = yUb;*/
		leftChild = ConstructSpatialTree(x,y,xLb,xLb+countXLeft,yLb,yLb+countYLeft, c, d, el, 1, node); 
	}
	else {
		assert(0); //should not come here.
	}
#endif
/*	node = new Vertex();
	node->parent = parent;
	node->level = parent->level + 1;
	leftChild = ConstructSpatialTree(x,y,xLb,xLb+countXLeft,yLb,yLb+countYLeft, c, d, el, 1, node); */

	c = center + rad/2; //if(!flag) then the value of c.
	if(flag) {
		c = (mid + (center + rad))/2; //overwrite the value of c if flag is true.
		d = ((center + rad) - mid)/2; 
	}
	
	assert(node != NULL);
	rightChild = ConstructSpatialTree(x,y,xLb+countXLeft,xUb,yLb+countYLeft,yUb, c, d, 0, er, node); 
	node->left = leftChild;
	node->right = rightChild;
	node->center = center;
	node->radius = rad;
	//node->eta = eta0/rad;
	return node;
}

/*This is a 1D problem. So, when the input points are sorted and when the interval that they fall into is bisected,  the resulting two sub-intervals can be identified using the index of the left-most and right-most point falling in the respective sub-interval. In addition to this trivial case, two nodes are classified as neighbors if their resp. center's are within a certain threshold distance.*/
bool AreAdjacent(const Vertex* node1, const Vertex* node2) {
	bool ret = true;
	assert((node1!=NULL) && (node2!=NULL));
		//trivial case.
		//if((node1->xRight == node2->xLeft) || (node1->xLeft == node2->xRight) || (node1->yRight == node2->yLeft) || (node1->yLeft == node2->yRight))
		//is-within-a-threshold case.
		if((node1->radius + node2->radius) <= (TAU * std::abs(node1->center - node2->center)))
			ret=false;
	return ret;
}

void UpdateNeighbors(Vertex* node) {
	
	Vertex* parent = node->parent;

	//Add siblings as neighbors
	if(parent->left && (parent->left != node))
		node->neighbors.push_back(parent->left);
	
	if(parent->right && (parent->right != node))
		node->neighbors.push_back(parent->right);
#ifdef DEBUG
	if(parent->left && (parent->left != node))
		node->nbrs.push_back(parent->left->label);
	
	if(parent->right && (parent->right != node))
		node->nbrs.push_back(parent->right->label);
#endif

	//Add cousins as neighbors (if the center of cells represented by cousins are within a threshold distance)
	for(int i=0;i<parent->neighbors.size();i++) {	
		Vertex* parNeighbor = parent->neighbors[i];
		//For adaptive trees, parent's neighbor i.e. uncle may not have any children i.e. the current node, 'node', may not have any cousins. In such a scenario, there is a possibility of uncle becoming a neighbor.
		if(!parNeighbor->isLeaf) {
			//check if any of the cousins are neighbors.
			if(AreAdjacent(parNeighbor->left, node))
				node->neighbors.push_back(parNeighbor->left);
			else {
				node->wellSeparatedNodes.push_back(parNeighbor->left);
				parNeighbor->left->wellSeparatedNodes.push_back(node);
			}
			
			if(AreAdjacent(parNeighbor->right, node))
				node->neighbors.push_back(parNeighbor->right);
			else {
				node->wellSeparatedNodes.push_back(parNeighbor->right);
				parNeighbor->right->wellSeparatedNodes.push_back(node);
			}
#ifdef DEBUG
			if(AreAdjacent(parNeighbor->left, node))
				node->nbrs.push_back(parNeighbor->left->label);
			else{
				node->il.push_back(parNeighbor->left->label);
				parNeighbor->left->il.push_back(node->label);
			}

			if(AreAdjacent(parNeighbor->right, node))
				node->nbrs.push_back(parNeighbor->right->label);
			else {
				node->il.push_back(parNeighbor->right->label);
				parNeighbor->right->il.push_back(node->label);
			}
#endif
		} 
		else {
			//check if uncle is a neighbor. If so, add it to node's neighbor list
			if(AreAdjacent(node, parNeighbor)) {
				node->neighbors.push_back(parNeighbor);
				if(node->isLeaf)
					parNeighbor->neighbors.push_back(node);
			}
			else {
				node->wellSeparatedNodes.push_back(parNeighbor);
				parNeighbor->wellSeparatedNodes.push_back(node);
			}
#ifdef DEBUG
			if(AreAdjacent(parNeighbor, node)) {
				node->nbrs.push_back(parNeighbor->label);
				if(node->isLeaf)
					parNeighbor->nbrs.push_back(node->label);
			} else {
				node->il.push_back(parNeighbor->label);
				parNeighbor->il.push_back(node->label);
			}
#endif

		}
	}
	

	if(node->left)
		UpdateNeighbors(node->left);
	if(node->right)
		UpdateNeighbors(node->right);
	return;
}


void PreorderVisitor(Vertex* node, TreeVisitor* vis) {
		vis->PreorderActions(node);
		if(node->left)
			PreorderVisitor(node->left, vis);

		if(node->right)
			PreorderVisitor(node->right, vis);

	return;
}



void PostorderVisitor(Vertex* node, TreeVisitor* vis)
{
		if(node->left)
			PostorderVisitor(node->left, vis);

		if(node->right)
			PostorderVisitor(node->right, vis);

		vis->PostorderActions(node);
	return;
}

void PostorderVisitorLeaf(Vertex* node, TreeVisitor* vis)
{
		if(node->left)
			PostorderVisitorLeaf(node->left, vis);

		if(node->right)
			PostorderVisitorLeaf(node->right, vis);

		vis->LeafNodeActions(node);
	return;
}

/* Goal: to accelerate a matrix vector multiplication Phi*q. Phi is a 'separable' matrix that can be created using column vectors x and y. 
*fun is the kernel (e.g. used to determine electrostatic potential, gravitational force etc.)
*r is the number of terms used in the Taylor series expansion. 
* org specifies the order of elements in x for determining gap.
* gap is the difference between consecutive elements in a global order. E.g. if we have interlacing vectors x and y i.e. x[i]<y[i]<x[i+1]<y[i+1], gap[i]=y[i]-x[i]
* fun is the kernel used. E.g. 1/(x-y) or 1/(x-y)^2 etc.
* By default, scaling is on (scaling is used for stability) */
double* fmm1d_local_shift_2(int r, double *x, double *y, double * q, const double *p_gap, const int* p_org, const int fun, const int numXElems, int numYElems, const int scaling=1) {
	/*dataX=x;
	dataY=y;
	dataQ=q;
	scalingLocal=scaling;
	funLocal=fun;
	numTerms=r;
	org=p_org;
	gap=p_gap;
	orgSize=numXElems;*/
	std::qsort(x, numXElems, sizeof(double), compare);  
	std::qsort(y, numYElems, sizeof(double), compare); 
	double z2 = std::max(x[numXElems-1], y[numYElems-1]);
	double z1 = std::min(x[0], y[0]);
	z2 += 0.1 * std::abs(z2);
	z1 -= 0.1 * std::abs(z1);

	Vertex* markerNode = new Vertex(); // this is a dummy node created to avoid having 'if(parent) == NULL' checks in ConstructSPatialTree.
	markerNode->level=0;
	Vertex* node = ConstructSpatialTree(x,y,0,numXElems,0,numYElems, (z1+z2)/2, (z2-z1)/2, 1, 1, markerNode); 
	markerNode->left=node;
#ifdef DEBUG
	AssignLabels(node);
#endif

	if(node->left)
		UpdateNeighbors(node->left);

	if(node->right)
		UpdateNeighbors(node->right);

	TreeVisitor* v=new TreeVisitor(x, y, q, r, scaling, fun, p_org, p_gap, numXElems, node);
	PostorderVisitor(node,v);
	PreorderVisitor(node, v);
	v->result=new double[numXElems];
	for(int i=0;i<numXElems;i++) v->result[i]=0;
	PostorderVisitorLeaf(node, v);
	double* ret=v->result;
	delete v; //v->result is not destroyed.
	delete markerNode;
	return ret;
}	

void TriFMM1TreeVisitor::PreorderActions(Vertex* node) {
	if(node->level >= 2) {
		for(std::vector<Vertex*>::iterator iter=node->wellSeparatedNodes.begin();iter!=node->wellSeparatedNodes.end();iter++) {
			double* outB=NULL; 
			Vertex* wsNode=*iter;
			//printf("ComputeB called with eta0:%lf a:%lf b:%lf, dx:%lf, dy:%lf, fun:%d, scaling:%d\n",eta0, node->center, wsNode->center, 2*(node->radius),2*(wsNode->radius),funLocal, scalingLocal); 
			ComputeB_Scaled(&outB, numTerms, eta0, node->center, wsNode->center, 2*(node->radius), 2*(wsNode->radius), funLocal, scalingLocal);
			if(node->center > wsNode->center) {
				assert(wsNode->res);
				//lower triangular part
				if(node->computed == false) {
					//computing ul{i}=Bij * v{j}
					node->ul=new double[numTerms];
					for(int j=0;j<numTerms;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += (outB[j*numTerms+i] * wsNode->res[i]); 
						}
						node->ul[j] = temp;
					}
					node->computed = true;
				} else {
					//computing ul{i}+=Bij * v{j}
					for(int j=0;j<numTerms;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += (outB[j*numTerms+i] * wsNode->res[i]); 
						}
						node->ul[j] += temp;
					}
				}
			} else if(node->center < wsNode->center) {
				//upper triangular part
				if(node->computedUpper == false) {
					//computing uu{i}=Bij * v{j}
					node->uu = new double[numTerms];
					for(int j=0;j<numTerms;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += (outB[j*numTerms+i] * wsNode->res[i]); 
						}
						node->uu[j] = temp;
					}
					node->computedUpper = true;
				} else {
					//computing uu{i}+=Bij * v{j}
					for(int j=0;j<numTerms;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += (outB[j*numTerms+i] * wsNode->res[i]); 
						}
						node->uu[j] += temp;
					}
				}
			}
			else {assert(0);}
			delete [] outB;
		}
	

		if(node->level > 3) {
		    Vertex* parent = node->parent;
		    if(parent->computed || parent->computedUpper) {
			double *outT=NULL;
			ComputeT_Scaled(&outT, numTerms, eta0, node->center, parent->center,  2*(node->radius), 2*(parent->radius), scalingLocal);
			if(parent->computed) {
				assert(parent->ul);
				if(node->computed == false) {
					//computing ul{i}=Ri * ul{p}
					node->ul=new double[numTerms];
					for(int j=0;j<numTerms;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += (outT[j*numTerms+i] * parent->ul[i]); 
						}
						node->ul[j] = temp;
					}
					node->computed = true;
				} else {
					//computing ul{i}+=Ri * ul{p}
					for(int j=0;j<numTerms;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += (outT[j*numTerms+i] * parent->ul[i]); 
						}
						node->ul[j] += temp;
					}
				}
			}
			
			if(parent->computedUpper) {
				assert(parent->uu);
				if(node->computedUpper == false) {
					//computing uu{i}=Ri * uu{p}
					node->uu=new double[numTerms];
					for(int j=0;j<numTerms;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += (outT[j*numTerms+i] * parent->uu[i]); 
						}
						node->uu[j] = temp;
					}
					node->computedUpper=true;
				} else {
					//computing u{i}+=Ri * u{p}
					for(int j=0;j<numTerms;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += (outT[j*numTerms+i] * parent->uu[i]); 
						}
						node->uu[j] += temp;
					}
				}
			}
			delete [] outT;
		    }
		}	    
	}
}

void TriFMM1TreeVisitor::LeafNodeActions(Vertex* node) {
	if((node->left==NULL) && (node->right==NULL)) {
		if((node->xRight - node->xLeft) != 0) {
			double* outU=NULL;

			int numElements=node->xRight-node->xLeft;
			if(numElements > 0) {
				ComputeU_Scaled(&outU, dataX+node->xLeft, numElements, numTerms, eta0, node->center, 2*(node->radius), scalingLocal);
				GetTransposeInPlace(outU, numTerms, numElements);
				if(node->computed) {
					//compute zl(px, :) = Ui * ul{i};
					assert(node->zl == NULL);
					node->zl=new double[numElements];
					for(int j=0;j<numElements;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += outU[j*numTerms+i] * node->ul[i]; 
						}
						node->zl[j] = temp; 
						assert((node->xLeft+j) < orgSize);
						result[node->xLeft+j] = temp;
					}
				}
				if(node->computedUpper) {
					//compute zu(px, :) = Ui * uu{i};
					assert(node->zu == NULL);
					node->zu=new double[numElements];
					for(int j=0;j<numElements;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += outU[j*numTerms+i] * node->uu[i]; 
						}
						node->zu[j] = temp; 
						assert((orgSize+node->xLeft+j) < 2*orgSize);
						result[orgSize+node->xLeft+j] = temp;
					}
				}
				delete [] outU;
			}
			else {
				/*if(node->computed) {
						delete [] node->res;
						node->res = new double[orgSize];
						for(int i=0;i<numElements;i++)
							node->res[i]=0;
				}
				if(node->computedUpper) {
						delete [] node->resu;
						node->resu = new double[orgSize];
						for(int i=0;i<numElements;i++)
							node->resu[i]=0;
				}*/
				assert(0);
			}
			
			
			int xiVectorSize = node->xRight-node->xLeft;
			assert(xiVectorSize);
			double* ti = new double[xiVectorSize];
			double* yi = new double[xiVectorSize];
			for(int i=0;i<xiVectorSize;i++) {
				assert((node->xLeft+i) < orgSize);
				assert(i < orgSize);
				yi[i] = dataY[org[node->xLeft+i]]; 
				ti[i] = gap[node->xLeft+i];
			}

			//Union(i, neighbor{i})
			bool isSelf=true; 
			node->neighbors.insert(node->neighbors.begin(), node);
			for(std::vector<Vertex*>::iterator iter=node->neighbors.begin();iter!=node->neighbors.end();iter++) {
				Vertex* nbr = *iter;
				int xjVectorSize = nbr->yRight-nbr->yLeft;
				if(xjVectorSize !=0) {
					double *D = new double[xiVectorSize*xjVectorSize];
					//computing D=yi-yj.' 
					for(int i=0;i<xiVectorSize;i++) {
						for(int j=0;j<xjVectorSize;j++) {
					 		D[i*xjVectorSize+j] = yi[i] - dataY[nbr->yLeft+j];
						}
					}

					bsxfun('p',&D, std::make_pair(xiVectorSize, xjVectorSize), ti, std::make_pair(xiVectorSize, 1)); 
					switch(funLocal) {
						case 1: 
							{
								for(int i=0;i<xiVectorSize;i++) 
									for(int j=0;j<xjVectorSize;j++) 
					 					D[i*xjVectorSize+j] = 1/(double)D[i*xjVectorSize+j];
							}
							break;
						case 2:
							{
								for(int i=0;i<xiVectorSize;i++) 
									for(int j=0;j<xjVectorSize;j++) 
					 					D[i*xjVectorSize+j] = 1/(double)(D[i*xjVectorSize+j] * D[i*xjVectorSize+j]);
							}
							break;
						case 3:
							{
								for(int i=0;i<xiVectorSize;i++) { 
									for(int j=0;j<xjVectorSize;j++) {
					 					D[i*xjVectorSize+j] = std::log(std::abs(D[i*xjVectorSize+j]));
										if(isSelf && isinf(D[i*xjVectorSize+j]))
											D[i*xjVectorSize+j] = 0;
									}
								}
							}
							break;
						default:
							assert(0);				
					}
					
					double *Dl = new double[xiVectorSize*xjVectorSize];
					memcpy(Dl, D, sizeof(double)*xiVectorSize*xjVectorSize);
					const double* px=dataX+node->xLeft;
					const double* py=dataY+nbr->yLeft;
					for(int i=0;i<xiVectorSize;i++) {
						for(int j=0;j<xjVectorSize;j++) {
							if(px[i] < py[j] )
					 			Dl[i*xjVectorSize+j] = 0 ;
						}
					}
					
					for(int i=0;i<xiVectorSize;i++)
						for(int j=0;j<xjVectorSize;j++)
							D[i*xjVectorSize+j]=D[i*xjVectorSize+j]-Dl[i*xjVectorSize+j];

					//zl(px, :) = zl(px, :) + Dl * q(py, :);
                    			//zu(px, :) = zu(px, :) + D * q(py, :);
					//cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, xiVectorSize, 1, xjVectorSize, 1, Dl, xjVectorSize, dataQ+nbr->yLeft, 1, 1, node->res, 1);
					//cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, xiVectorSize, 1, xjVectorSize, 1, D, xjVectorSize, dataQ+nbr->yLeft, 1, 1, node->resu, 1);
					cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, xiVectorSize, 1, xjVectorSize, 1, Dl, xjVectorSize, dataQ+nbr->yLeft, 1, 1, result+node->xLeft, 1);
					cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, xiVectorSize, 1, xjVectorSize, 1, D, xjVectorSize, dataQ+nbr->yLeft, 1, 1, result+orgSize+node->xLeft, 1);	
					/*//memcpy(result+node->xLeft, node->res, sizeof(double)*xiVectorSize);
					//memcpy(result+orgSize+node->xLeft, node->resu, sizeof(double)*xiVectorSize);
					for(int i=0;i<xiVectorSize;i++) 
						result[node->xLeft+i] = node->res[i];
					for(int i=0;i<xiVectorSize;i++) 
						result[orgSize+node->xLeft+i] = node->resu[i];*/
					delete [] D;
					delete [] Dl;
				}
				isSelf=false; //set flag to false in all iterations except the first one. When this is false, (isnan() call is not done)
			}
			delete [] yi;
			delete [] ti;
			node->neighbors.erase(node->neighbors.begin());
		}
	}
	
}

double* trifmm1d_local_shift(int r, double *x, double *y, double * q, const double *p_gap, const int* p_org, const int fun, const int numXElems, int numYElems, const int scaling=1) {
	/*dataX=x;
	dataY=y;
	dataQ=q;
	scalingLocal=scaling;
	funLocal=fun;
	numTerms=r;
	org=p_org;
	gap=p_gap;
	orgSize=numXElems;*/
	std::qsort(x, numXElems, sizeof(double), compare);  
	std::qsort(y, numYElems, sizeof(double), compare); 
	double z2 = std::max(x[numXElems-1], y[numYElems-1]);
	double z1 = std::min(x[0], y[0]);
	z2 += 0.1 * std::abs(z2);
	z1 -= 0.1 * std::abs(z1);

	Vertex* markerNode = new Vertex(); // this is a dummy node created to avoid having 'if(parent) == NULL' checks in ConstructSPatialTree.
	markerNode->level=0;

	Vertex* node = ConstructSpatialTree(x,y,0,numXElems,0,numYElems, (z1+z2)/2, (z2-z1)/2, 1, 1, markerNode); 
	markerNode->left=node;
#ifdef DEBUG
	AssignLabels(node);
#endif
	/*std::ofstream of;//("trifmmoutput_afterpreorder.txt",std::ofstream::out);
	PrintTree(node, of);*/
	

	if(node->left)
		UpdateNeighbors(node->left);

	if(node->right)
		UpdateNeighbors(node->right);

	TreeVisitor* v=new TriFMM1TreeVisitor(x, y, q, r, scaling, fun, p_org, p_gap, numXElems, node);
	PostorderVisitor(node,v);
	PreorderVisitor(node, v);
	//PrintTree(node, of);
        v->result = new double[numXElems*2];
	for(int i=0;i<numXElems*2;i++) v->result[i]=0;
	PostorderVisitorLeaf(node, v); 
	double* ret=v->result;
	delete v; // v->result is not destroyed
	delete markerNode;
	return ret;
}

void FMM1TreeVisitor::LeafNodeActions(Vertex* node) {
	if((node->left==NULL) && (node->right==NULL)) {
		if((node->xRight - node->xLeft) != 0) {
			double* outU=NULL;
			int numElements=node->xRight-node->xLeft;
			if(numElements > 0) {
			/*PrintArray<int>(&numElements, 1, 1, "fmminput_Uscaled.txt");
			PrintArray<const double>(dataX+node->xLeft, 1, numElements, "fmminput_Uscaled.txt");
			PrintArray<int>(&numTerms, 1, 1, "fmminput_Uscaled.txt");
			PrintArray<double>(&eta0, 1, 1, "fmminput_Uscaled.txt");
			double arg1=node->center;
			double arg2=2*(node->radius);
			PrintArray<double>(&arg1, 1, 1, "fmminput_Uscaled.txt");
			PrintArray<double>(&arg2, 1, 1, "fmminput_Uscaled.txt");
			PrintArray<int>(&scalingLocal, 1, 1, "fmminput_Uscaled.txt");*/


				ComputeU_Scaled(&outU, dataX+node->xLeft, numElements, numTerms, eta0, node->center, 2*(node->radius), scalingLocal);
				//PrintArray2<double>(outU, numTerms, numElements, "fmmoutput_Uscaled_BT.txt");
				GetTransposeInPlace(outU, numTerms, numElements);
				//PrintArray2<double>(outU, numElements, numTerms, "fmmoutput_Uscaled.txt");
				//assert(0);
				if(node->computed) {
					//compute z(px, :) = Ui * u{i};
					node->zl=new double[numElements];
					for(int j=0;j<numElements;j++){
						double temp=0;
						for(int i=0;i<numTerms;i++){
							temp += outU[j*numTerms+i] * node->ul[i]; 
						}
						node->zl[j] = temp;
						assert((node->xLeft+j) < orgSize);
					       result[node->xLeft+j] = temp; 	
					}
				}
				delete [] outU;
			} else{
				assert(0);
				/*if(node->computed) {
						delete [] node->res;
						node->res = new double[orgSize];
						for(int i=0;i<numElements;i++)
							node->res[i]=0;
				}*/
			}
			int xiVectorSize = node->xRight-node->xLeft;
			assert(xiVectorSize);
			double* ti = new double[xiVectorSize];
			double* yi = new double[xiVectorSize];
			for(int i=0;i<xiVectorSize;i++) {
				assert((node->xLeft+i) < orgSize);
				assert(i < orgSize);
				yi[i] = dataY[org[node->xLeft+i]]; 
				ti[i] = gap[node->xLeft+i];
			}


			//Union(i, neighbor{i})
			bool isSelf=true; 
			node->neighbors.insert(node->neighbors.begin(), node);
			for(std::vector<Vertex*>::iterator iter=node->neighbors.begin();iter!=node->neighbors.end();iter++) {
				Vertex* nbr = *iter;
				int xjVectorSize = nbr->yRight-nbr->yLeft;
				if(xjVectorSize !=0) {
					double *D = new double[xiVectorSize*xjVectorSize];
					//computing D=yi-yj.' 
					for(int i=0;i<xiVectorSize;i++) {
						for(int j=0;j<xjVectorSize;j++) {
							assert((nbr->yLeft+j) < orgSize);
					 		D[i*xjVectorSize+j] = yi[i] - dataY[nbr->yLeft+j];
						}
					}
					bsxfun('p',&D, std::make_pair(xiVectorSize, xjVectorSize), ti, std::make_pair(xiVectorSize, 1)); 
					switch(funLocal) {
						case 1: 
							{
								for(int i=0;i<xiVectorSize;i++) 
									for(int j=0;j<xjVectorSize;j++) 
					 					D[i*xjVectorSize+j] = 1/(double)D[i*xjVectorSize+j];
							}
							break;
						case 2:
							{
								for(int i=0;i<xiVectorSize;i++) 
									for(int j=0;j<xjVectorSize;j++) 
					 					D[i*xjVectorSize+j] = 1/(double)(D[i*xjVectorSize+j] * D[i*xjVectorSize+j]);
							}
							break;
						case 3:
							{
								for(int i=0;i<xiVectorSize;i++) { 
									for(int j=0;j<xjVectorSize;j++) {
					 					D[i*xjVectorSize+j] = std::log(std::abs(D[i*xjVectorSize+j]));
										if(isSelf && isinf(D[i*xjVectorSize+j]))
											D[i*xjVectorSize+j] = 0;
									}
								}
							}
							break;
						default:
							assert(0);				
					}

					//cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, xiVectorSize, 1, xjVectorSize, 1, D, xjVectorSize, dataQ+nbr->yLeft, 1, 1, node->res, 1);
					cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, xiVectorSize, 1, xjVectorSize, 1, D, xjVectorSize, dataQ+nbr->yLeft, 1, 1, result+node->xLeft, 1);
                    			//z(px, :) = z(px, :) + D * q(py, :);
					delete [] D;
				}
				isSelf=false; //set flag to false in all iterations except the first one. When this is false, (isnan() call is not done)
			}

			delete [] yi;
			delete [] ti;
			node->neighbors.erase(node->neighbors.begin());
		}
	}
	
}


double* fmm1d_local_shift(int r, double *x, double *y, double * q, const double *p_gap, const int* p_org, const int fun, const int numXElems, int numYElems, const int scaling=1) {
	/*dataX=x;
	dataY=y;
	dataQ=q;
	scalingLocal=scaling;
	funLocal=fun;
	numTerms=r;
	org=p_org;
	gap=p_gap;
	orgSize=numXElems;*/
	std::qsort(x, numXElems, sizeof(double), compare);  
	std::qsort(y, numYElems, sizeof(double), compare); 
	double z2 = max(x[numXElems-1], y[numYElems-1]);
	double z1 = min(x[0], y[0]);
	z2 += 0.1 * abs(z2);
	z1 -= 0.1 * abs(z1);

	Vertex* markerNode = new Vertex(); // this is a dummy node created to avoid having 'if(parent) == NULL' checks in ConstructSPatialTree.
	markerNode->level=0;

	Vertex* node = ConstructSpatialTree(x,y,0,numXElems,0,numYElems, (z1+z2)/2, (z2-z1)/2, 1, 1, markerNode); 
	markerNode->left=node;
#ifdef DEBUG
	AssignLabels(node);
#endif
	

	if(node->left)
		UpdateNeighbors(node->left);

	if(node->right)
		UpdateNeighbors(node->right);

	TreeVisitor* v=new FMM1TreeVisitor(x, y, q, r, scaling, fun, p_org, p_gap, numXElems, node);
	PostorderVisitor(node,v);
	PreorderVisitor(node, v);
	/*std::ofstream of;//("fmmoutputtree_u.txt",std::ofstream::out);
	PrintTree(node, of);*/
        v->result = new double[numXElems];
	for(int i=0;i<numXElems;i++) v->result[i]=0;
	PostorderVisitorLeaf(node, v);
	double* ret=v->result;
	delete v; //v->result is not destroyed
	delete markerNode;
	return ret;
}	


#ifdef DEBUG
void AssignLabels(Vertex* node) {
	int label=1;
	//assign labels as per level ordered sequence of traversing the tree.
	std::queue<Vertex*> Q;
	Q.push(node);
	while(!Q.empty()) {
		Vertex* node=Q.front();	
		Q.pop();
		node->label = label++;
		if(node->left)
			Q.push(node->left);
		if(node->right)
			Q.push(node->right);
	}
}
#endif
