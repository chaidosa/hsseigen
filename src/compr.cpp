#include<string.h>
#include<algorithm>
#include<cmath>
#include<stdio.h>
#include<assert.h>
#include"compr.h"
#include"QR.h"

/*A is the input parameter initialized to T[i] for leaf i.
Q is the output parameter computing Ui, R is the output parameter computing Ri */
void compr(double* A, std::pair<int, int>aSize, double** Q, std::pair<int, int>& qSize, double** R, std::pair<int, int>& rSize, char* tol, double par)
{
	int opt = 4;
	if((aSize.first == 0) || (aSize.second == 0))
		return;

	switch(opt)
	{
		case 4:
			{
				//allocate space for Q. = numRows * numRows;
				*Q = new double[aSize.first*aSize.first];
				qSize.first  = aSize.first;
				qSize.second = aSize.first;
				*R = new double[aSize.first*aSize.second];
				rSize  = aSize;
				//allocate space for pivot = numCols;
				int *P = new int[aSize.second];
				memset(P,0,sizeof(int)*aSize.second);
				qpr(*Q,*R,P,A, aSize.first, aSize.second);	
				int rk=-1; 
				if(!strcmp(tol, "tol"))
				{
					/*get the last index for which abs(diag(R)/R[1,1]) >= par  */
					double r11 = **R;
					for(int i=0;i<rSize.first;i++)
					{
						double diagR = *((*R)+(i*rSize.second)+i); 
						double tmp = diagR/r11;
						if(tmp < 0)
							tmp *= -1;
						if(tmp > par)
							rk = i;		
					}
					rk++;
				}
				else
					rk = par;
				
				rk = std::min(rk,rSize.second);
				rk = std::min(rk,rSize.first);
				
				//truncate Q and R. keep only rk rows in Q and R. Keep all the original clumns in Q and only the pivots in R.
				//The pivots are returned in P at the end of qpr.
				double* tmpQ = new double[rk*qSize.first];
				for(int i=0;i<qSize.first;i++)
					memcpy(tmpQ+(i*rk), (*Q)+(i*qSize.second),sizeof(double)*rk);
				qSize.second=rk;
	
				delete [] *Q;
				*Q=tmpQ;
				
#ifdef DEBUG
				for(int i=0;i<aSize.second;i++)
				{
					if((P[i] < 1) || (P[i]>aSize.second))
					{
						printf("Invalid column ID\n");
						assert(0);
					}
				}
#endif
				//cannot think of an efficient way to copy columns referred to by the elements of P vector.
				int* RP = new int[aSize.second];
				for(int i=0;i<aSize.second;i++)
					RP[P[i]-1]=i;

				rSize.first = rk;
				double *tmpR = new double[rk*rSize.second];
				for(int i=0;i<rk;i++)
					for(int j=0;j<rSize.second;j++)
						tmpR[(i*rSize.second)+j]=(*R)[i*rSize.second+RP[j]];
				delete [] RP;
				delete [] *R;
				*R=tmpR;
				delete [] P;
			}
			break;
		default:
			printf("Compression not supported yet\n");
			break;
	}
}
