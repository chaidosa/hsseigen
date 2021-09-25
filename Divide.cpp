#include<string.h>
#include<assert.h>
#include "cblas.h"
#include "QR.h"
#include "Divide.h"


DivideOutParams* Divide(HSSMat* hssMat, BinTree* tr)
{
	DivideOutParams* ret=new DivideOutParams();
	int n2 = tr->numNodes;
	const int *lv = tr->level;   
	int r = 0;     //r will be HSS block-rank of matrix, initialized to 0
	std::vector<int> cl = tr->leaves; //array of leaf level nodes
	
	//this array holds the children ids of every node in the tree.
	std::vector<std::pair<int, int> > ch(n2);
	//index is array of leaf level nodes within list of all nodes
	int *index =new int[n2] ;
	memset(index, 0, sizeof(int)*n2);
	for(int i=0;i<n2;i++)     
    	{
		std::vector<int> children = tr->GetChildren(i+1);
		if(children.empty())
        	{
			index[i]=i+1;
			ch[i].first = -1;
			ch[i].second = -1; 
		}
		else
		{
			ch[i].first = children[0];
			ch[i].second = children[1]; 
		}
    	}

	int n3 = cl.size();
	int* m = new int[n3];
	int n=0;	//size of matrix
	for(int i=0; i<n3;i++)
	{
		m[i] = hssMat->uSizes[cl[i]-1].first;
		n += m[i];
	}

	std::vector<int> l(n3+1);
	l[0]=1;
	for(int i=0;i<n3;i++)
		l[i+1]=l[i]+m[i];

	int nb =hssMat->N -1; //Check for B's that are NULL.
	for(int i=0;i<nb;i++)
	{
		if(hssMat->B[i])
		{
			std::pair<int,int> dim = hssMat->bSizes[i];
			r = std::max(r,std::max(dim.first, dim.second));
		}
	}

	//array of siblings
	int* sib = new int[n2];
	memset(sib,0,sizeof(int)*n2);
	for(int i=0;i<n2;i++)
	    sib[i] = tr->GetSibling(i+1);

	int mlv=tr->maxLevel;
	
	//allocate memory for rank-1 updates
	double *Z= new double[n*r*mlv];
	memset(Z,0,sizeof(double)*n*r*mlv);
	std::pair<int, int> zSize(n,r*mlv);
	
	//determine coordinates of each HSS block. 
	std::vector<std::pair<int, int> > lvec(n2);
	int z=0;
	for(int i=0;i<n2;i++)
	{
		std::pair<int, int> children = ch[i];
		if(children.first==-1)
		{
			z = z+1;
			lvec[i].first = l[z-1];
			lvec[i].second = l[z]-1;
		}
    		else
		{
			//Not a leaf. Get the leftmost and rightmost leafNodes of the subtree.
			int aa  = tr->GetLeftMostChild(i+1);
			int bb  = tr->GetRightMostChild(i+1);
			lvec[i].first = lvec[aa-1].first;
			lvec[i].second = lvec[bb-1].second;
		}
	}

	// divide. traverse tree, top-down
	int d = mlv;
	for(int k = 0;k<mlv;k++)
	{
		d = d - 1;
		//visit each node at level k
		for(int ck =0; ck<n2;ck++) 
		{
			if((lv[ck] == k-1) && (ch[ck].first != -1))
		    	{
	
				int ca = ch[ck].first; // first child of node ck
		    		int cb = ch[ck].second;// second child of node ck
		    		std::vector<int> des1 = tr->GetDescendents(ca); // descendents of first child of node ck
		    		std::vector<int> des2 = tr->GetDescendents(cb); // descendents of second child of node ck
			    	int s1 = des1.size();
			    	int s2 = des2.size();
			    	double* Tf1 = NULL;
			    	double* Tf2 = NULL;
				std::pair<int, int> tf1Size(0,0), tf2Size(0,0);

				int ii;
			    	for(int i=0;i<s2;i++)
				{
					//form rank-r updates
					int j=((s2-1) < i)?(s2-1):i;
					ii = des1[i];
					int jj = des2[j];
					int tmpm = (hssMat->uSizes[ii-1]).first;
					int tmpk = (hssMat->uSizes[ii-1]).second;
					int tmpn = (hssMat->rSizes[ii-1]).second;
					assert(tmpk == (hssMat->rSizes[ii-1]).first);
					double* T1 =new double[tmpm*tmpk];
					memcpy(T1, hssMat->U[ii-1],sizeof(double)*tmpm*tmpk);
					double* B = hssMat->R[ii-1];
					double* A = T1;
					double* C = new double[tmpm*tmpk];
					double alpha=1.0, beta=0;
					while(lv[ii-1] > lv[jj-1])
					{
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, tmpm, tmpn, tmpk, alpha, A, tmpk, B, tmpn, beta, C, tmpn);
						memcpy(A,C,sizeof(double)*tmpm*tmpk);
						ii = tr->tr[ii-1];
					}

					int _tmpm = (hssMat->uSizes[jj-1]).first;
					int _tmpk = (hssMat->uSizes[jj-1]).second;
					int _tmpn = (hssMat->rSizes[jj-1]).second;
					assert(_tmpk == (hssMat->rSizes[jj-1]).first);
					double* T2 =new double[_tmpm*_tmpk];
					memcpy(T2, hssMat->U[jj-1],sizeof(double)*_tmpm*_tmpk);
					double* D = T2;
					double* E = hssMat->R[jj-1];
					double* F = new double[_tmpm*_tmpk];
					while(tr->tr[ii-1] != tr->tr[jj-1])
					{
						//performs T1=T1*R{ii}
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, tmpm, tmpn, tmpk, alpha, A, tmpk, B, tmpn, beta, C, tmpn);
						memcpy(A,C,sizeof(double)*tmpm*tmpk);
						ii = tr->tr[ii-1];

						//performs T2=T2*R{jj}
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _tmpm, _tmpn, _tmpk, alpha, D, _tmpk, E, _tmpn, beta, F, _tmpn);
						memcpy(D,F,sizeof(double)*_tmpm*_tmpk);
				    		jj = tr->tr[jj-1];
					}
					
					delete [] C;
					delete [] F;

					//Tf1 is derived from T1
					double* newTf1 = new double[(tf1Size.first+tmpm)*(tf1Size.second)];
					if(Tf1)
					{
						assert(tf1Size.second == tmpk);
						memcpy(newTf1, Tf1, sizeof(double)*tf1Size.first*tf1Size.second);
					}
					int hRowIndex=0;
					for(hRowIndex=0;hRowIndex<tmpm;hRowIndex++)
						memcpy(newTf1+((hRowIndex+tf1Size.first)*tf1Size.second),T1+(tmpk*hRowIndex),sizeof(double)*tf1Size.second);
					tf1Size.first += tmpm;
					tf1Size.second = tmpk;
					if(Tf1)
						delete [] Tf1;
					Tf1=newTf1;
					
					//Tf2 is derived from T2
					if((s2-1) >= i)
					{
						//Tf1 is derived from T1
						double* newTf2 = new double[(tf2Size.first+_tmpm)*(tf2Size.second)];
						if(Tf2)
						{
							assert(tf2Size.second == _tmpk);
							memcpy(newTf2, Tf2, sizeof(double)*tf2Size.first*tf2Size.second);
						}
						int hRowIndex=0;
						for(hRowIndex=0;hRowIndex<_tmpm;hRowIndex++)
							memcpy(newTf2+((hRowIndex+tf2Size.first)*tf2Size.second),T2+(_tmpk*hRowIndex),sizeof(double)*tf2Size.second);
						tf2Size.first += _tmpm;
						tf2Size.second = _tmpk;
						if(Tf2)
							delete [] Tf2;
						Tf2=newTf2;
					}

					delete [] T1;
					delete [] T2;
				}
			    	// store rank-one updates
				int rt = tf1Size.second;
				//Z(l(1,index(des1(1,1),1)):l(1,index(des1(1,s1),1)+1)-1,d*r+1:d*r+rt) = Tf1;
				int rSI = l[index[des1[0]-1]-1]-1;
				int rEI = l[index[des1[s1-1]-1]]-2;
				int cSI = d*r;
				int cEI = d*r+rt-1;
				int numRows = rEI-rSI+1;
				int numCols = cEI-cSI+1;
				assert(numCols == tf1Size.second);
				for(int tmpRowIndex=0;tmpRowIndex<numRows;tmpRowIndex++)
					memcpy(Z+((tmpRowIndex+rSI)*zSize.second)+cSI,Tf1+(tmpRowIndex*tf1Size.second),sizeof(double)*numCols);


				//Z(l(1,index(des2(1,1),1)):l(1,index(des2(1,s2),1)+1)-1,d*r+1:d*r+rt) = Tf2*B{i1}';
				assert((hssMat->bSizes[ii-1]).second == tf2Size.second);
				//first get Tf2*B{ii}'
				int tmpBRowWidth = (hssMat->bSizes[ii-1]).second;
				int tmpBColWidth = (hssMat->bSizes[ii-1]).first;
				double* transMatrix = new double[tmpBRowWidth * tmpBColWidth];
				memcpy(transMatrix, hssMat->B[ii-1], sizeof(double)*tmpBRowWidth*tmpBColWidth);
				GetTransposeInPlace(transMatrix,tmpBColWidth, tmpBRowWidth);
				double* interMatrix = new double[tf2Size.first * tmpBColWidth];
				//performs Tf2*B{ii}'
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, tf2Size.first, tmpBColWidth, tf2Size.second, 1.0, Tf2, tf2Size.second, transMatrix, tmpBColWidth, 0.0, interMatrix, tmpBColWidth);
				delete [] transMatrix;
				//now copy interMatrix to the mentioned submatrix of Z.
				rSI = l[index[des2[0]-1]-1]-1;
				rEI = l[index[des2[s2-1]-1]]-2;
				cSI = d*r;
				cEI = d*r+rt-1;
				numRows = rEI-rSI+1;
				numCols = cEI-cSI+1;
				assert(numCols == tmpBColWidth);
				for(int tmpRowIndex=0;tmpRowIndex<numRows;tmpRowIndex++)
					memcpy(Z+((tmpRowIndex+rSI)*zSize.second)+cSI,interMatrix+(tmpRowIndex*tmpBColWidth),sizeof(double)*numCols);
				delete [] interMatrix;


				delete [] Tf1;
				delete [] Tf2;
			}
		}

		/*double** D = new double*[tr->leaves.size()];
		memset(D,0,sizeof(double*)tr->leaves.size()];
		std::vector<std::pair<int, int> >& dSizes=ret->dSize;
		dSize.reserve(tr->leaves.size());
		for(int i=0;i<tr->leaves.size();i++)
			dSize[i]=std::make_pair(0,0);*/
		
		double** D=hssMat->D;
		

		// update D generators
		for(int j=0; j<n2;j++)
		{
			std::pair<int, int> curDSize = hssMat->dSizes[j];
        		if(ch[j].first == -1)
			{
           			//D[j] = D[j] - Z(l(1,index(j,1)):l(1,index(j,1)+1)-1,d*r+1:d*r+r)*Z(l(1,index(j,1)):l(1,index(j,1)+1)-1,d*r+1:d*r+r)';
			
				//first get part of Z (tmpZ), and tmpZ'
				int rSI = l[index[j]-1]-1;
				int rEI = l[index[j]]-2;
				int cSI = d*r;
				int cEI = d*r+r-1;
				int numRows = rEI-rSI+1;
				int numCols = cEI-cSI+1;
				double* tmpZ = new double[numRows*numCols];

				//Z
				assert(numRows == numCols);
				assert(numRows == curDSize.first);
				assert(numCols == curDSize.second);
				for(int tmpRowIndex=0;tmpRowIndex<numRows;tmpRowIndex++)
					memcpy(tmpZ+tmpRowIndex*numCols,Z+((tmpRowIndex+rSI)*zSize.second)+cSI,sizeof(double)*numCols);
				//Z'
				double* tmpZTrans = new double[numRows*numCols]; 
				memcpy(tmpZTrans, tmpZ, sizeof(double)*numRows*numCols);
				GetTransposeInPlace(tmpZTrans,numRows, numCols);

				//dgemm computes C=alpha*op(A)*op(B)+beta*C. set alpha=-1. beta=1 A=Z B=Z' C=D m=numRows 
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numRows, numCols, numCols, -1.0, tmpZ, numCols, tmpZTrans, numRows, 1.0, D[j], numRows);
        		}

        		//update B generators
			int j2 = tr->GetSibling(j+1);
			if(lv[j] > k)
			{
			    if((j2==-1)||(j < j2-1))
			    {
				int an = j+1;
				while(lv[an-1] > k)
				    an = tr->tr[an-1];

				int i1, i2, sib_an=tr->GetSibling(an);
				if(an < sib_an)
				{
				   i1 = an;
				   i2 = sib_an;
				}
				else
				{
				    i2 = an;
				    i1 = sib_an;
				}

				double* B_new1 = hssMat->R[j];
				double* tmpBnew1 = new double[hssMat->rSizes[j].first*hssMat->rSizes[j].second];
		
				int numRowsBnew2=hssMat->rSizes[j2-1].first;
				int numColsBnew2=hssMat->rSizes[j2-1].second;
				double* B_new2 = new double[numRowsBnew2*numColsBnew2];
				memcpy(B_new2,hssMat->R[j2-1],sizeof(double)*numRowsBnew2*numColsBnew2);
				GetTransposeInPlace(B_new2,numRowsBnew2, numColsBnew2);
				std::swap(numRowsBnew2, numColsBnew2);
				double* tmpBnew2 = new double[numRowsBnew2*numColsBnew2];

				double* transB = new double[hssMat->bSizes[i1-2].first*hssMat->bSizes[i1-2].second];
				memcpy(transB,hssMat->B[i1-2],sizeof(double)*hssMat->bSizes[i1-2].first*hssMat->bSizes[i1-2].second);
				GetTransposeInPlace(transB,hssMat->bSizes[i1-2].first, hssMat->bSizes[i1-2].second);

				int a = j+1;
				if(a < i2) 
				{
				    a = tr->tr[a-1];
				    int _a = j+1;
				    int _m,_n,_k;
				    while(lv[_a-1] - lv[i1-1] > 1)
				    {
					    _a = tr->tr[_a-1];
					    double* transRA = new double[hssMat->rSizes[_a-1].first*hssMat->rSizes[_a-1].second];
					    memcpy(transRA, hssMat->R[_a-1],hssMat->rSizes[_a-1].first*hssMat->rSizes[_a-1].second);
					    GetTransposeInPlace(transRA,hssMat->rSizes[_a-1].first, hssMat->rSizes[_a-1].second);
					    int numRowsTransRA=hssMat->rSizes[_a-1].second, numColsTransRA=hssMat->rSizes[_a-1].first;

					    //B_new1 = B_new1 * hssMat->R[a-1];
					    //dgemm computes C=alpha*op(A)*op(B)+beta*C. set alpha=1. beta=0.0 A=B_new1 B=hssMat->R C=tmpBnew1 m=hssMat->rSizes[j].first, k=hssMat->rSizes[j].second, n=hssMat->rSizes[a-1].second, lda=k, ldb=n, ldc=n
					    _m=hssMat->rSizes[j].first;_n=hssMat->rSizes[a-1].second;_k=hssMat->rSizes[j].second;
					    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _m, _n, _k, 1.0, B_new1, _k, hssMat->R[a-1], _n, 0.0, tmpBnew1, _n);
					    B_new1 = tmpBnew1;
					    //B_new2 = transRA * B_new2;
					    //dgemm computes C=alpha*op(A)*op(B)+beta*C. set alpha=1. beta=0.0 A=transRA B=B_new2 C=tmpBnew2 m=hssMat->rSizes[j].first, k=hssMat->rSizes[j].second, n=hssMat->rSizes[a-1].second, lda=k, ldb=n, ldc=n
					    _m=numRowsTransRA;_n=numColsBnew2;_k=numColsTransRA;
					    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _m, _n, _k, 1.0, transRA, _k, B_new2, _n, 0.0, tmpBnew2, _n);
					    B_new2 = tmpBnew2;
					    delete [] transRA;
				    }

				    if(a == i2)
				    {
					//B_new1 = B_new1 * hssMat->B[i1-1]';
					_m=hssMat->rSizes[j].first;_n=hssMat->bSizes[i1-2].first;_k=hssMat->rSizes[j].second;
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _m, _n, _k, 1.0, B_new1, _k, transB, _n, 0.0, tmpBnew1, _n);
					B_new1 = tmpBnew1;

					//B_new2 = hssMat->B[i1-1] * B_new2;
					_m=hssMat->bSizes[i1-2].first;_n=numColsBnew2;_k=hssMat->bSizes[i1-2].second;
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _m, _n, _k, 1.0, hssMat->B[i1-2], _k, B_new2, _n, 0.0, tmpBnew2, _n);
					B_new2 = tmpBnew2;
				    }
				    //hssMat->B[j] = hssMat->B[j] - B_new1*B_new2;
				    //dgemm computes C=alpha*op(A)*op(B)+beta*C. set alpha=-1.0 beta=1.0 A=B_new1 B=B_new2 C=hssMat->B[j-1] m=hssMat->rSizes[j].first, k=hssMat->rSizes[j].second, n=numColsBnew2, lda=k, ldb=n, ldc=n
				    _m=hssMat->rSizes[j].first;_n=numColsBnew2;_k=numColsBnew2;
				    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _m, _n, _k, -1.0, B_new1, _k, B_new2, _n, 1.0, hssMat->B[j], _n);
				}
				
				delete [] transB;
				delete [] tmpBnew1;
				delete [] tmpBnew2;
				delete [] B_new1;
				delete [] B_new2;
			    }
			}
		}
	}

	ret->Z = Z;
	ret->zSize = zSize;
	ret->l = l;
	ret->lvec = lvec;
	delete [] index;
	delete [] m;
	delete [] sib;
}
