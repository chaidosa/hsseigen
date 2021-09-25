#include<string.h>
#include<math.h>
#include<algorithm>
#include<assert.h>
#include"mat2hsssym.h"
#include "compr.h"
#include "QR.h"

/* Computes Generators, D, U, R, B for nodes of the binary tree. 
 * Leaf nodes have only U and R matrices computed (by a call to compr; dependent on T[i], where i is the leaf node) */
HSSMat* mat2hsssym(double* A, int aSize, BinTree* bt, int* m, int mSize, char* tol, double par)  {
	HSSMat* ret = new HSSMat();
	int n = mSize; //n is number of leaves
	int n0=n;
	int lSize=0;
	int aColWidth, aRowWidth = sqrt(aSize); // since A is a symmetric matrix, aRowWidth = aCols and (aRowWidth*aCols = aSize)
	assert((aRowWidth * aRowWidth) == aSize);
	aColWidth = aRowWidth;
	int N = bt->GetNumNodes();
	assert(bt->GetLeaves().size() == (N+1)/2);
	int numLeaves = (N+1)/2;
	assert(numLeaves == n);

	int* l = new int[mSize+1];
	lSize = mSize+1;
	l[0]=1;
	for(int i=0;i<mSize;i++)
		l[i+1]=l[i]+m[i];
	
	//create a list to hold the matrix dimensions of the diagonal matrices 
	std::pair<int, int>* dSizes = new std::pair<int, int>[numLeaves];
	for(int i=0;i<bt->leaves.size();i++)
		dSizes[i]=std::make_pair(0,0);

	//create a temporary array to hold the entire row-block except the diagonal block of A.
	double** T=new double*[N];
	std::pair<int, int>* tSizes = new std::pair<int, int>[N];
	for(int i=0;i<N;i++)
		tSizes[i]=std::make_pair(0,0);

	//copy the diagonal block of A to D.
	double** D=new double*[numLeaves];
	for(int i=0;i<n0;i++)
	{
		int ii=bt->leaves[i];
		int dRowWidth = l[i+1]-l[i];
		int dColWidth = l[i+1]-l[i];
		assert(dRowWidth == dColWidth);
		D[i] = new double[dRowWidth * dColWidth];
		dSizes[i] = std::make_pair(dColWidth, dRowWidth);
		double * tmp = D[i];
		//row and column start and end indices.
		int rSI = l[i]-1;
		int rEI = l[i+1]-1;
		int cSI = rSI;
		int cEI = rEI;
		//copy the A submatrix in [rSI, cSI), [rEI, cEI) to D
		for(int k=0,j=rSI;j<rEI;j++,k++)
			memcpy(tmp+(k*dRowWidth), A+j*aRowWidth+cSI, sizeof(double)*dRowWidth);

		int tColWidth = dColWidth;
		int tRowWidth = aRowWidth-dRowWidth;
		T[ii-1] = new double[tRowWidth * tColWidth];
		double* tmpT = T[ii-1];
		tSizes[ii-1] = std::make_pair(tColWidth, tRowWidth);
		
		if(cSI == 0)
		{
			for(int k=0,j=rSI;j<rEI;j++,k++)
				memcpy(tmpT+(k*tRowWidth), A+j*aRowWidth+cEI, sizeof(double)*tRowWidth);
		}
		else if((cSI > 0) && (cEI < aRowWidth))
		{
			for(int k=0,j=rSI;j<rEI;j++,k++)
			{
				memcpy(tmpT+(k*tRowWidth), A+j*aRowWidth, sizeof(double)*(cSI));
				memcpy(tmpT+(k*tRowWidth)+cSI, A+j*aRowWidth+cEI, sizeof(double)*(aRowWidth-cEI));
			}
		}
		else if(cEI == aRowWidth)
		{
			for(int k=0,j=rSI;j<rEI;j++,k++)
				memcpy(tmpT+(k*tRowWidth), A+j*aRowWidth, sizeof(double)*tRowWidth);
		}
		else
		{
			printf("ERROR in creating T Matrix rSI:%d rEI:%d cSI:%d cEI:%d\n",rSI,rEI,cSI,cEI);
			assert(0);
		}
	}
	
	/*double *Q1 = new double[100*100];
	double *R1 = new double[100*3100];
	int *P1 = new int[3100];
	memset(P1,0,sizeof(int)*3100);
	qpr(Q1, R1, P1, T[0], 100, 3100);*/
	
	//create array bn.
	std::vector<int> bn(n);
	for(int i=0;i<n;i++)
		bn[i]=i;	

	//Ceate temporary array for bn.
	std::vector<int> bn0 = bn;


	//Ceate temporary arrays for l and m.
	int* l0 = new int[lSize];
	int* m0 = new int[mSize];
	int tmpMSize = mSize, tmpLSize = lSize;
	

	std::vector<int> cl, cl0 = bt->leaves;
	cl = cl0;
	//U and B
	//create space for one U and B matrix for every node of the tree except root.
	double** U=new double*[N-1];
	memset(U, 0, sizeof(double*)*(N-1));
	std::pair<int, int>* uSizes = new std::pair<int, int>[N-1];
	for(int i=0;i<N-1;i++)
		uSizes[i] = std::make_pair(0,0);
	
	double** B=new double*[N-1];
	memset(B, 0, sizeof(double*)*(N-1));
	std::pair<int, int>* bSizes = new std::pair<int, int>[N-1];
	for(int i=0;i<N-1;i++)
		bSizes[i] = std::make_pair(0,0);

	for(int i=0;i<N-1;i++)
	{
		//get the index of node being updated
		std::vector<int>::iterator iter = std::find(cl.begin(), cl.end(),i+1);
		assert(iter != cl.end());
		int ibn = iter - cl.begin();
		int ii = bn[ibn];

		memcpy(l0,l,sizeof(int)*tmpLSize);
		memcpy(m0,m,sizeof(int)*tmpMSize);
		
		std::vector<int> ch = bt->GetChildren(i+1);
		if(ch.size() ==0)
		{
			//leaf node.
			double* tmpT = T[i];
			std::pair<int, int> uSize(0,0);
			double* R=NULL;
			std::pair<int, int> rSize=uSize;
			printf("LeafNode: Computing U%d R%d\n",i+1,i+1); 
			compr(tmpT, tSizes[i], &(U[i]),uSize, &R, rSize, tol,par);
			delete [] tmpT;
			T[i] = R;
			tSizes[i]=rSize;
			uSizes[i]=uSize;
		}
		else
		{
			//get the left child index = c1
			std::vector<int>::iterator iter = std::find(cl0.begin(), cl0.end(),ch[0]);
			assert(iter != cl0.end());
			ibn = iter - cl0.begin();
			int i1 = bn0[ibn];
			//get the right child index = c2
			iter = std::find(cl0.begin(), cl0.end(),ch[1]);
			assert(iter != cl0.end());
			ibn = iter - cl0.begin();
			int i2 = bn0[ibn];

			double* H=NULL;
			std::pair<int, int> hSize;
			int ch0RowWidth = tSizes[ch[0]-1].second;	
			int ch0ColWidth = tSizes[ch[0]-1].first;	
			int ch1RowWidth = tSizes[ch[1]-1].second;	
			int ch1ColWidth = tSizes[ch[1]-1].first;	
			
			double* ch0Matrix = T[ch[0]-1];
			double* ch1Matrix = T[ch[1]-1];
			if(i1==1)
			{
				//Case 1: H = T[ch[0]] columns l(i1+2)-m(i1) to end. +  T[ch[1]] columns l(i1+2)-m(i2) to end
				int cSI = l[i1+2]-m[i1]-1;			
				int hColWidth = ch0ColWidth + ch1ColWidth;
				int hRowWidth = ch0RowWidth - cSI; 
				assert(hRowWidth == (ch1RowWidth-l[i1+2]+m[i2]+1));
				hSize=std::make_pair(hColWidth, hRowWidth);
				H = new double[hSize.first*hSize.second];

				int hRowIndex=0;
				for(hRowIndex=0;hRowIndex<ch0ColWidth;hRowIndex++)
					memcpy(H+(hRowIndex*hRowWidth),ch0Matrix+(ch0RowWidth*hRowIndex)+cSI,sizeof(double)*hRowWidth);
			
				cSI = l[i1+2]-m[i2]-1;
				for(int ch1RowIndex=0;ch1RowIndex<ch1ColWidth;hRowIndex++, ch1RowIndex++)
					memcpy(H+(hRowIndex*hRowWidth),ch1Matrix+(ch1RowWidth*ch1RowIndex)+cSI,sizeof(double)*hRowWidth);
				
			}
			else if(i2==bn.size()-1)
			{
				//Case 2: H = T[ch[0]] columns 1 to l(i1)-1 +  T[ch[1]] columns 1 to l(i1)-1
				int cSI = 0;	
				int hColWidth = ch0ColWidth + ch1ColWidth;
				int hRowWidth = l[i1]-1; //1 to l[i1]-1 => l[i1] columns
				hSize=std::make_pair(hColWidth, hRowWidth);
				H = new double[hSize.first*hSize.second];

				int hRowIndex=0;
				for(hRowIndex=0;hRowIndex<ch0ColWidth;hRowIndex++)
					memcpy(H+(hRowIndex*hRowWidth),ch0Matrix+(ch0RowWidth*hRowIndex)+cSI,sizeof(double)*hRowWidth);
			
				for(int ch1RowIndex=0;ch1RowIndex<ch1ColWidth;hRowIndex++, ch1RowIndex++)
					memcpy(H+(hRowIndex*hRowWidth),ch1Matrix+(ch1RowWidth*ch1RowIndex)+cSI,sizeof(double)*hRowWidth);

			}
			else
			{
				//Case 3: combination of cases 1 and 2. H = T[ch[0]] columns 1 to l(i1)-1 and  l(i1+2)-m(i1) to end. +  
				//T[ch[1]] columns 1 to l(i1)-1 and l(i1+2)-m(i2) to end.
				int cSI = l[i1+2]-m[i1]-1;			
				int hColWidth = ch0ColWidth + ch1ColWidth;
				int hRowWidth = ch0RowWidth - cSI + l[i1]-1; 

				hSize=std::make_pair(hColWidth, hRowWidth);
				H = new double[hSize.first*hSize.second];

				int hRowIndex=0;
				int remainingRowWidth = ch0RowWidth - cSI;
				for(hRowIndex=0;hRowIndex<ch0ColWidth;hRowIndex++)
				{
					cSI=0;
					memcpy(H+(hRowIndex*hRowWidth),ch0Matrix+(ch0RowWidth*hRowIndex)+cSI,sizeof(double)*(l[i1]-1));
					cSI = l[i1+2]-m[i1]-1;
					memcpy(H+(hRowIndex*hRowWidth)+(l[i1]-1),ch0Matrix+(ch0RowWidth*hRowIndex)+cSI,sizeof(double)*remainingRowWidth);
				}
			
				assert(remainingRowWidth == ch1RowWidth-l[i1+2]+m[i2]+1);
				for(int ch1RowIndex=0;ch1RowIndex<ch1ColWidth;hRowIndex++, ch1RowIndex++)
				{
					cSI=0;
					memcpy(H+(hRowIndex*hRowWidth),ch1Matrix+(ch1RowWidth*ch1RowIndex)+cSI,sizeof(double)*(l[i1]-1));
					cSI = l[i1+2]-m[i2]-1;
					memcpy(H+(hRowIndex*hRowWidth)+(l[i1]-1),ch1Matrix+(ch1RowWidth*ch1RowIndex)+cSI,sizeof(double)*remainingRowWidth);
				}
			}
			
	

			//B Matrix
			int bColWidth = ch0ColWidth;
			int bRowWidth = (l[i1+2]-m[i1]-1) - (l[i1+1]-m[i1]) + 1;
			bSizes[ch[0]-1]=std::make_pair(bColWidth, bRowWidth);
			
			printf("Computing B%d\n",ch[0]-1); 
			B[ch[0]-1] = new double[bColWidth * bRowWidth];
			double *tmpBMatrix = B[ch[0]-1];
			int cSI = (l[i1+1]-m[i1]-1);
			for(int ch0RowIndex=0;ch0RowIndex<ch0ColWidth;ch0RowIndex++)
				memcpy(tmpBMatrix+(bRowWidth*ch0RowIndex),ch0Matrix+(ch0RowIndex*ch0RowWidth)+cSI,sizeof(double)*bRowWidth); 
			
			//U Matrix	
			std::pair<int, int> uSize(0,0);
			double* R=NULL;
			std::pair<int, int> rSize=uSize;
			printf("Computing U%d\n",i); 
			compr(H, hSize, &(U[i]), uSize, &R, rSize, tol, par);
			delete [] H;
			T[i] = R;
			tSizes[i] = rSize;
			uSizes[i] = uSize;

        		m[ii] = tSizes[i].first;
        		for(int tmpmIndex = ii+1;tmpmIndex < tmpMSize; tmpmIndex++)
				m[tmpmIndex] = m[tmpmIndex+1];
			tmpMSize--;

			//update partition intervals
			l[0]=1;
			for(int j=0;j<tmpMSize;j++)
				l[j+1]=l[j]+m[j];
			tmpLSize--;
		}
		
		/*for every j in cl - contains all leaves initially and then internal nodes as leaves get updated.
		{
			if(j != current matrix being updated, ii)
			{
				//j1 = get the ID of the node whose matrix is being updated =cl(j)
				//k1 = level of current-level matrix tree node
				//Get j1's T Matrix and initialze H. H = T{j1}
				//if(the node corresponding to current matrix being updated has leftsiblings, uncles, granduncles, and so on) //if(ii>1)
				//  for every j > ii (left sibling, uncle, granduncle etc.)  //if( ii > j)
				//    Intialize T{j1} with either 0 or from H (1 to l(ii)-m(j)-1) columns.
				//  else
				//    Initialize T{j1} from H (1 to l(ii)-1) columns.
				//else
				//  initialize T{j1}=0;
				//if(j is granduncle, uncle, left sibling)
				//  update T{j1} with columns from T{i} (l(j) to l(j+1)-1 columns transposed), where i is the current node being updated (top-level)
				//else
				//  update T{j1} with columns from T{i} (l(j)-m(ii) to l(j+1)-m(ii)-1 columns transposed), where i is the current node being updated (top-level)
				//
				//if current matrix being updated, ii,  is < n, where n is 32
				//  if(node corersponding to current top-level matrix is leaf //empty(ch)
				//    if j < ii (left sibling, uncle, granduncle)
				//      Update T{j1} with columns from H (l(ii+1)-m(j):end)
				//    else
				//      Update T{j1} with columns from H (l(ii+1):end)
				//  else
				//    if j < ii (left sibling, uncle, granduncle)
				//      Update T{j1} with columns from H (l0(ii+2)-m0(j):end)
				//    else
				//      Update T{j1} with columns from H (l0(ii+2):end)
				//    
				//      
			}	
		}*/

		/* if(top-level matrix node is leaf ) //empty(ch)
 			update partition size, m(ii) = size(T(i),1)
		        update entire l array.	   	
 		*/ 
	
		/* prevcl = currentcl cl0=cl
 		   if(i is rightchild)
			remove left and right children from cl and add their parent's ID
			prevbn = currentbn
			update currentbn from rightchildID to end (shift left the elements and reduce by 1)
 		*/ 	

#if 1
		for(int j=0;j<cl.size();j++)
		{
			if(j != ii)
			{
				int j1 = cl[j];
				//printf("%d:%d ",j+1,j1);
				double* H = T[j1-1];
				std::pair<int, int> hSize = tSizes[j1-1];
				if(ii > 0)
				{
					if(ii > j)
					{
						if((j==0) && (ii==1))
						{
							T[j1-1]=NULL;	
							tSizes[j1-1] = std::make_pair(0, 0);
						}
						else
						{
							//T[j1-1]=H(:,1:l(ii)-m(j)-1)
							int tmpRowWidth = l[ii]-m[j]-1;
							double *tmp = new double[hSize.first*tmpRowWidth];
							for(int tmpRowIndex=0;tmpRowIndex<hSize.first;tmpRowIndex++)
								memcpy(tmp+(tmpRowIndex*tmpRowWidth),H+(tmpRowIndex*hSize.second), sizeof(double)*tmpRowWidth);
							T[j1-1]=tmp;
							tSizes[j1-1] = std::make_pair(hSize.first, tmpRowWidth);
						}
					}
					else
					{
						//T[j1-1]=H(:,1:l(ii)-1)
						int tmpRowWidth = l[ii]-1;
						double *tmp = new double[hSize.first*tmpRowWidth];
						for(int tmpRowIndex=0;tmpRowIndex<hSize.first;tmpRowIndex++)
							memcpy(tmp+(tmpRowIndex*tmpRowWidth),H+(tmpRowIndex*hSize.second), sizeof(double)*tmpRowWidth);
						T[j1-1]=tmp;
						tSizes[j1-1] = std::make_pair(hSize.first, tmpRowWidth);
					}	
				}
				else
				{
					T[j1-1]=NULL;
					tSizes[j1-1] = std::make_pair(0, 0);
				}

				double* transMatrix=NULL;
				int transMatrixRowWidth=0, transMatrixColWidth=0;
				int cSI, tMatrixRowWidth = tSizes[i].second;
				transMatrixColWidth=tSizes[i].first;
				if(j < ii)
				{
					//T[j1-1] = [T[j1-1], T[i](:,l(j):l(j+1)-1)']
					transMatrixRowWidth=(l[j+1]-l[j]);
					cSI = l[j]-1;
				}
				else
				{
					//T[j1-1] = [T[j1-1], T[i](:,l(j)-m(ii):l(j+1)-m(ii)-1)']
					transMatrixRowWidth=((l[j+1]-m[ii]) - (l[j]-m[ii]));
					cSI = l[j]-m[ii]-1;
				}

				transMatrix = new double[transMatrixRowWidth*transMatrixColWidth];
				for(int tmpRowIndex=0;tmpRowIndex<transMatrixColWidth;tmpRowIndex++)
					memcpy(transMatrix+(tmpRowIndex*transMatrixRowWidth),T[i]+(tmpRowIndex*tMatrixRowWidth)+cSI,sizeof(double)*transMatrixRowWidth);
				//now get the transpose of transMatrix.
				GetTransposeInPlace(transMatrix, transMatrixColWidth, transMatrixRowWidth);
				std::swap(transMatrixRowWidth, transMatrixColWidth);



				//now append all columns of transMatrix to T[j1-1].
				double* tmp1=NULL, *tmp = T[j1-1];
				int tmpColWidth=transMatrixColWidth, tmpRowWidth=0, tmp1RowWidth=0, remainingRowWidth;
				if(tmp)
				{
					tmp1=tmp;
					tmp1RowWidth = tSizes[j1-1].second;
				}
				cSI = 0;
				remainingRowWidth = transMatrixRowWidth - cSI;
				tmpRowWidth = tSizes[j1-1].second + remainingRowWidth;
				tmp = new double[tmpRowWidth* tmpColWidth];
				for(int tmpRowIndex =0;tmpRowIndex<tmpColWidth;tmpRowIndex++)
				{
					if(tmp1)
						memcpy(tmp+(tmpRowIndex*tmpRowWidth),tmp1+(tmpRowIndex*tmp1RowWidth), sizeof(double)*tmp1RowWidth);
					memcpy(tmp+(tmpRowIndex*tmpRowWidth)+tmp1RowWidth,transMatrix+(tmpRowIndex*transMatrixRowWidth)+cSI, sizeof(double)*remainingRowWidth);
				}

				delete [] transMatrix;
				if(tmp1)
					delete [] tmp1;
				T[j1-1]=tmp;
				tSizes[j1-1] = std::make_pair(tmpColWidth, tmpRowWidth);


				//Append specific columns of HMatrix to T[j1-1].
				if(ii < n)
				{
					tmp1=NULL;
					tmpColWidth=hSize.first;
					int hRowWidth = hSize.second;
					if(tmp)
					{
						tmp1=tmp;
						tmp1RowWidth = tSizes[j1-1].second;
					}

					if(ch.size() == 0)
					{
						if(j < ii) 		//T{j1-1}= [T{j1-1} H(:,l(ii+1)-m(j):end]
							cSI = (l[ii+1]-m[j]) - 1;
						else 			//T{j1-1}= [T{j1-1} H(:,l(ii+1)):end]
							cSI = l[ii+1] - 1;
					}
					else
					{
						if(j < ii)		//T{j1-1}= [T{j1-1} H(:,l0(ii+2)-m0(j):end]
							cSI = (l0[ii+2]-m0[j]) - 1;
						else			//T{j1-1}= [T{j1-1} H(:,l0(ii+2)):end]
							cSI = l0[ii+2] - 1;
					}

					remainingRowWidth = hRowWidth - cSI;
					tmpRowWidth = tSizes[j1-1].second + remainingRowWidth;
					tmp = new double[tmpColWidth* tmpRowWidth];
					for(int tmpRowIndex =0;tmpRowIndex<tmpColWidth;tmpRowIndex++)
					{
						if(tmp1)
							memcpy(tmp+(tmpRowIndex*tmpRowWidth),tmp1+(tmpRowIndex*tmp1RowWidth), sizeof(double)*tmp1RowWidth);
						memcpy(tmp+(tmpRowIndex*tmpRowWidth)+tmp1RowWidth,H+(tmpRowIndex*hRowWidth)+cSI, sizeof(double)*remainingRowWidth);
					}

					if(tmp1)
						delete [] tmp1;
					T[j1-1]=tmp;
					tSizes[j1-1] = std::make_pair(tmpColWidth, tmpRowWidth);
				}
				delete [] H;
			}
		}

#endif
		//printf("%d \n",ii+1);
		/*for(int t=0;t<N;t++)
			printf("%d (%d %d),  ",t,tSizes[t].first,tSizes[t].second);
		printf("\n");*/
		if(ch.size() == 0)
		{
			m[ii] = tSizes[i].first;

			//update partition intervals
			l[0]=1;
			for(int j=0;j<tmpMSize;j++)
				l[j+1]=l[j]+m[j];
		}

		cl0=cl;
		ch = bt->GetChildren(bt->tr[i]);
		if(ch[1] == i+1)
		{
			std::vector<int>::iterator iter = std::find(cl0.begin(),cl0.end(),ch[0]);
			assert(iter!=cl0.end());
			int ibn1 = iter-cl0.begin();
			iter = std::find(cl0.begin(),cl0.end(),ch[1]);
			assert(iter!=cl0.end());
			int ibn2 = iter-cl0.begin();
			cl[bn0[ibn1]] = bt->tr[i];
			cl.erase(cl.begin()+bn0[ibn2]);
			
			bn0=bn;
			bn.erase(bn.begin()+ibn2);
			for(int j=ibn2;j<bn.size();j++)
				bn[j] -=1;
		}
	}

	//B matrix of 1st level.
	int i = N, i1, ibn;
	std::vector<int> ch = bt->GetChildren(i);
	std::vector<int>::iterator iter = std::find(cl0.begin(),cl0.end(),ch[0]);
	assert(iter!=cl0.end());
	ibn = iter-cl0.begin(); i1 = bn0[ibn];

	//B{ch(1)} = T{ch(1)}(:,l(i1+1)-m(i1):l(i1+2)-m(i1)-1);
	int bRowWidth = (l[i1+2]-m[i1])- (l[i1+1]-m[i1]);
	int bColWidth = tSizes[ch[0]-1].first;
	bSizes[ch[0]-1]=std::make_pair(bColWidth,bRowWidth);
	int ch0RowWidth = tSizes[ch[0]-1].second;
	assert(B[ch[0]-1] == NULL);
	printf("Updating B%d\n",ch[0]-1); 
	B[ch[0]-1] = new double[bRowWidth * bColWidth];
	double *tmpBMatrix = B[ch[0]-1], *ch0Matrix=T[ch[0]-1];
	int cSI = (l[i1+1]-m[i1]-1);
	for(int ch0RowIndex=0;ch0RowIndex<bColWidth;ch0RowIndex++)
		memcpy(tmpBMatrix+(ch0RowIndex * bRowWidth),ch0Matrix+(ch0RowIndex*ch0RowWidth)+cSI,sizeof(double)*bRowWidth); 


	//R matrix
	double** R = new double*[N-1];
	memset(R,0,sizeof(double*)*(N-1));
	std::pair<int, int>* rSizes = new std::pair<int, int>[N-1];
	for(int i=0;i<N-1;i++)
		rSizes[i] = std::make_pair(0,0);

	for(int i=N-2;i>=0;i--)
    	{
		ch = bt->GetChildren(i+1);
    		if(ch.size() > 0)
		{
			int sz = uSizes[ch[0]-1].second;
			//R[ch[0]-1] = U{i}(1:sz,:);
			int rMatrixColWidth = sz, rMatrixRowWidth = uSizes[i].second;
			printf("Updating R%d\n",ch[0]-1); 
			R[ch[0]-1] = new double[rMatrixColWidth * rMatrixRowWidth];
			double* tmp = R[ch[0]-1];
			for(int tmpRowIndex=0;tmpRowIndex<rMatrixColWidth;tmpRowIndex++)
				memcpy(tmp+(tmpRowIndex*rMatrixRowWidth), U[i]+(tmpRowIndex*rMatrixRowWidth),sizeof(double)*rMatrixRowWidth);
			rSizes[ch[0]-1] = std::make_pair(rMatrixColWidth,rMatrixRowWidth);
			//R[ch[1]-1] = U{i}(sz+1:end,:);
			rMatrixColWidth = uSizes[i].first-sz, rMatrixRowWidth = uSizes[i].second;
			printf("Updating R%d\n",ch[1]-1); 
			R[ch[1]-1] = new double[rMatrixColWidth * rMatrixRowWidth];
			tmp = R[ch[1]-1];
			for(int tmpRowIndex=0;tmpRowIndex<rMatrixColWidth;tmpRowIndex++)
				memcpy(tmp+(tmpRowIndex*rMatrixRowWidth), U[i]+((tmpRowIndex+sz)*rMatrixRowWidth),sizeof(double)*rMatrixRowWidth);
			rSizes[ch[1]-1] = std::make_pair(rMatrixColWidth,rMatrixRowWidth);
		}
	}

	for(int i=0;i<N;i++)
	{
		ch = bt->GetChildren(i+1);   
		if(ch.size() > 0) 
		{
			if(U[i])
				delete [] U[i];
			U[i] = NULL;
			uSizes[i]=std::make_pair(0,0);
		}
	}

	ret->D=D;ret->U=U;ret->R=R;ret->B=B;
	ret->dSizes=dSizes;ret->uSizes=uSizes;ret->rSizes=rSizes;ret->bSizes=bSizes;
	ret->N=N;

	return ret;
}
