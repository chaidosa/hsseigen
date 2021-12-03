//Calculating QR using SVD
#include<string.h>
#include<algorithm>
#include<cmath>
#include<stdio.h>
#include<assert.h>
#include"compr.h"
#include"QR.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
    #include<cblas.h>
}

/*A is the input parameter initialized to T[i] for leaf i.
Q is the output parameter computing Ui, R is the output parameter computing Ri */
void compr_new(double* A, std::pair<int, int>aSize, double** Q, std::pair<int, int>& qSize, double** R, std::pair<int, int>& rSize, char const* tol, double par)
{
    int minMN = std::min(aSize.first,aSize.second);
    double *S    = new double[minMN];
    double *U    = new double[aSize.first*(aSize.first)];
    double *VT   = new double[aSize.second*(aSize.second)];
    double *W	 = new double[minMN-1];

    //Note: m = aSize.first, n = aSize.second, lda = aSize.second, ldu = aSize.first, ldvt = aSize.second, info;
    //calling lapack routine for SVD singular value decomposition
    int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'A','A',aSize.first,aSize.second,A,aSize.second,S,U,aSize.first,VT,aSize.second,W);    
   if(info > 0) {
   	printf("ERROR: LAPACKE_dgesvd failed to converge\n");
	exit(1);
   }
    delete [] W;
    
    //To know the size of U, read On Return for u (ldu by atleast m): https://www.ibm.com/docs/en/essl/6.2?topic=llss-sgesvd-dgesvd-cgesvd-zgesvd-sgesdd-dgesdd-cgesdd-zgesdd-singular-value-decomposition-general-matrix
    *Q=U;
    qSize = {aSize.first,aSize.first}; 

    //find the index i in S where S[i] > par*S[0])
    double thresh=par*S[0];
    int rk=0; 
    for(rk=1;rk<minMN;rk++) {
	if(S[rk] <= thresh)
		break;
    }
    if(rk != aSize.first) {
    	//pick only the first rk columns of U. U is of size aSize.first x aSize.first (ROW_MAJOR).
	//When you look at U (now aSize.first x rk) the columns represent left singular vectors.
    	double *tempU = new double[aSize.first*rk];
	for(int i=0;i<aSize.first;i++)
		memcpy(tempU+i*rk,U+i*aSize.first,sizeof(double)*rk);
	*Q=tempU;
	qSize = {aSize.first,rk}; 
	delete [] U;

    	//pick only the first rk rows of VT. VT is of size aSize.second x aSize.second (ROW_MAJOR)
	//When you look at VT (now rk x aSize.second) the rows represent right singular vectors.
	double * tempVT = NULL;
	if( rk != aSize.second) {
		tempVT = new double[rk*aSize.second];
		for(int i=0;i<rk;i++)
			memcpy(tempVT+i*aSize.second,VT+i*aSize.second,sizeof(double)*aSize.second);
    		delete [] VT;
	}
	else
		tempVT = VT;
	
	//pick only the first rk singular values from S. Multiply S (representing a diagonal matrix of size rk x rk ) and tempVT ( matrix of size rk x aSize.second)
	for(int i=0;i<rk;i++)
		for(int j=0;j<aSize.second;j++)
			tempVT[i*aSize.second+j]=S[i]*tempVT[i*aSize.second+j];
	*R=tempVT;
    	rSize = {rk,aSize.second};
    }
    else {
	    //Output parameter Q is already set earlier. Compute output parameter R.
	    //R is computed as S*VT (S represents a diagonal matrix of size aSize.first x aSize.second of which only the min(aSize.first, aSize.second) diagonal elements are actually non-zeros. VT represents a matrix of size (aSize.second x aSize.second)
	    double *tempR = new double[aSize.first * aSize.second];
	    memset(tempR, 0, sizeof(double)*aSize.first*aSize.second);
	    if(aSize.first <= aSize.second) {
		    //rectangular matrix with number of rows smaller than number of cols. 
		    for(int i=0;i<aSize.first;i++)
			for(int j=0;j<aSize.second;j++)
				tempR[i*aSize.second+j]=S[i]*VT[i*aSize.second+j];
	    }
	    else {
		    //rectangular matrix with number of cols smaller than number rows. 
		    for(int i=0;(i<aSize.first) && (i<aSize.second);i++)
			for(int j=0;j<aSize.second;j++)
				tempR[i*aSize.second+j]=S[j]*VT[i*aSize.second+j];
	    }
	    *R = tempR;
    	    rSize = {aSize.first,aSize.second};
	    delete [] VT;
    }

    delete [] S;
}
