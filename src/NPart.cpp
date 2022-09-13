#include"NPart.h"
#include"BinTree.h"
#include<string.h>
#include<math.h>

/*Prerequisite (if number of nodes in the tree not specified i.e. ltr==0): n/ni indicates the number of leaves. Either this must be one less
 * than the maximum possible leaves in a complete binary tree or this must be exactly equal to the 
 * maximum number of leaves possible (a perfect power of two).
 * E.g. 23x23 matrix partitioned as 3 rows per partition i.e. n=23, ni=3. We would have 8 partitions with all the 
 * partitions except the last partition getting 3 rows of the matrix. Last partition gets 2 rows.
 * E.g. 44x44 matrix partitioned as 6 rows per partition i.e. n=44, ni=6, We would have 7 partitions with all the  
 * partitions except the last partition getting 6 rows of the matrix. Last partition gets 8 rows. */

void NPart(int n, int ni, BinTree** tr, int**m, int& mSize, int ltr)
{
	if((ltr==0) && (ni >= n)) {
		cout<<"Error partition size specified is greater than the number of rows in the input matrix!"<<endl;
		*tr=NULL;
		*m=NULL;
		return;
	}

	if(ltr == 0){
		//k = number of horizontal partitions of input matrix
		int k=n/ni;
		
		//create a list that contains the size of each partition
		*m  = new int[k];
		for(int i=0;i<k;i++)
			(*m)[i]=ni;

		/* After uniform partitioning, if the number of remainining rows of the input matrix contains more than half the size of a partition,
		 * create a new partition and assign the remainder of the rows to that partition */ 
		int mr = n%ni;
		if(mr >= ni/2)
		{
			k = k+1;	
			delete [] *m;
			*m  = new int[k];
			for(int i=0;i<k-1;i++)
				(*m)[i]=ni;
			(*m)[k-1]=mr;
		}
		else
			(*m)[k-1] += mr; //otherwise, add the remaining number of rows to the last partition already created.
		
		*tr =new BinTree(2*k-1);	
		if((*tr)->GetNumNodes() == 0)
		{
			//deleting the tree if not a complete binary tree.
			delete *tr;
			*tr = NULL;
		}

		//set number of leaf nodes (mSize)
		mSize = k;
	}
	else {
		/*if the number  of leaves are already specified, compute how many rows are to be assigned to each partition.
		 * number of partitions is equal to number of leaf nodes */
		int k=(ltr+1)/2;
		*m  = new int[k];
		ni = floor(n/k);
		//create and initialize a list that holds the sizes of all partitions
		for(int i=0;i<k;i++)
			(*m)[i]=ni;
		(*m)[k-1] += n % ni; //add the remaining number of rows to the last partition already created.
		mSize = k;
	}
}
