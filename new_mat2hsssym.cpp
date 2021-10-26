#include<string.h>
#include<math.h>
#include<stack>
#include<algorithm>
#include<assert.h>
#include"test_mat2hsssym.h"
#include "compr.h"
#include "QR.h"

using namespace std;
//assuming the 16x16 matrix
tHSSMat* t_mat2hsssym(double* A, int aSize, BinTree* bt, int* m, int mSize, char* tol, double par){
    tHSSMat* ret = new tHSSMat();
    int n = mSize; //number of leaf nodes
    int N = bt->GetNumNodes();
    int aRowWidth = sqrt(aSize);
    int aColWidth = aRowWidth; 
    
    //storing info about node starting and ending block  
    std::pair<int,int>* l = new std::pair<int,int>[N];
    for(int i=0;i<N;i++)
		l[i]=std::make_pair(0,0);    

    l[0]={0,m[0]-1};
    int it=0,lt=0;
    for(int i=0;i<N;i++)
    {
        std::vector<int> ch = bt->GetChildren(i+1);
        if(ch.size() == 0)
        {
            l[i] = {lt,lt+m[it]-1};
            lt = l[i].second +1;
            it=it+1;
        }
        else
        {
            l[i]={l[ch[0]-1].first,l[ch[1]-1].second};

        }
    }
    //To store the diagonal blocks
    double** D = new double*[N];

    //create a list to hold the matrix dimensions of the diagonal matrices 
	std::pair<int, int>* dSizes = new std::pair<int, int>[N];
    for(int i=0;i<N;i++)
		dSizes[i]=std::make_pair(0,0); 

    //to store generator U
    double** U = new double*[N];

    //Pair to hold the matrix dimension of generator U
    std::pair<int, int>* uSizes = new std::pair<int, int>[N];
    for(int i=0;i<N;i++)
		dSizes[i]=std::make_pair(0,0); 

    double** R = new double*[N];
    double** B = new double*[N];

    std::stack<vector<double>>S;
    //creating temmporary array toa hold entire row-block except the diagonal block
    double** T = new double*[N];
   
   //creating temporary list to hold dimension of temporary block
    std::pair<int, int>* tSizes = new std::pair<int, int>[N];
    for(int i=0;i<N;i++)
		tSizes[i]=std::make_pair(0,0); 

    //temp stack
    std::vector<int>St;
    int ns=0;
    //visiting every nodes in postorder 
    for(int i=0;i<N;i++){
        std::vector<int> ch = bt->GetChildren(i+1);
        //if the current node is a leaf node
        if(ch.size() == 0)
        {   
            int dRowWidth = l[i].second-l[i].first+1;
            int dColWidth = dRowWidth;
            D[i] = new double[dRowWidth*dColWidth];
            double *temp = D[i];
            dSizes[i]={dColWidth,dRowWidth};

            //row and column start indices
            int rSI = l[i].first;            
            int rEI = l[i].second;
            int cSI = rSI;
            int cEI = rEI;

            //creating or copying diagonal blocks from original matrix
            for(int k=0,j=rSI;j<=rEI;j++,k++)
                memcpy(temp+(k*dRowWidth),A+j*aRowWidth+cSI, sizeof(double)*dRowWidth);
            
            //dimension of row and column of temporary array 
            int tRowWidht = aColWidth-dColWidth;
            int tColWidht = dRowWidth;
            
            tSizes[i] = {tColWidht,tRowWidht}; 
            //size of the temporary matrix           
            T[i] = new double[tColWidht*tRowWidht];
            
            //temporary pointer to hold matrix data corresponding to index i.
            double* tempT = T[i];

            if (i==0)
            {
                for(int k=0,j=rSI;j<=rEI;j++,k++)
                    memcpy(tempT+(k*tRowWidht),A+j*aRowWidth+cEI+1, sizeof(double)*tRowWidht);

                std::pair<int,int>uSize(0,0);
                double* R = NULL;
                std::pair<int,int>rSize =uSize;
                printf("LeafNode: Computing U%d",i+1);
                compr(tempT, tSizes[i], &(U[i]),uSize, &R, rSize, tol, par);
                delete [] tempT;
                T[i]= R;
                tSizes[i]=rSize;
                uSizes[i]=uSize;
            }

            else
            {
                St.push_back(i-1);
                ns++;
                int current_pos;//used to keep track of the columns filled.
                for(int k=0;k<ns;k++)
                {
                    rSI = 0;
                    rEI = tSizes[St[k]].first-1;
                    cSI = tSizes[St[k]].second-(N-l[i].first)-1;
                    cEI = tSizes[St[k]].second-(N-l[i].second)-1;
                    double *tempTT = T[St[k]];
                    double* temp_array = new double[(rEI+1)*(cEI)]; 
                    int index=0;
                    for(int i=rSI;i<=rEI;i++){
                        for(int j=cSI;j<=cEI;j++){
                            temp_array[index++]=tempTT[j+i*tRowWidht];
                        }
                    }
                    
                    GetTransposeInPlace(temp_array,(rEI+1),(cEI+1));
                    for(int i=rSI,id=0;i<=rEI,id<index;i++){
                        for(int j=cSI;j<=cEI;j++){
                            tempT[j+i*tRowWidht]=temp_array[id++];
                            current_pos=j;
                        }
                    }
                    delete [] temp_array;                    
                }

                //copying remaining element to tempT array from A
                int cp = current_pos+1;
                rSI = l[i].first;
                rEI = l[i].second;
                cSI = l[i].second+1;
                cEI = aRowWidth;
                for(int i=rSI,tp=0;i<=rEI;i++,tp++){
                    for(int j=cSI;j<cEI;j++){
                        tempT[cp+tp*tRowWidht]=A[j+i*aRowWidth];
                        cp++;   
                    }
                    cp=current_pos+1;
                }

                std::pair<int,int>uSize(0,0);
                double* R = NULL;
                std::pair<int,int>rSize =uSize;
                printf("LeafNode: Computing U%d",i+1);
                compr(tempT, tSizes[i], &(U[i]),uSize, &R, rSize, tol, par);
                delete [] tempT;
                T[i]= R;
                tSizes[i]=rSize;
                uSizes[i]=uSize;

            }   
        }
        //if the currrent node is non-leaf node
        else
        {
                St.pop_back();
                ns--;
                int left = ch[0]-1;
                int right = ch[1]-1;
                int tColWidth = l[i].second-l[i].first +1;
                int tRowWidth = aRowWidth-tColWidth;
                T[i] = new double[tRowWidth*tColWidth];
                double* tempT = T[i];

            
             int r1 = tSizes[left].first;
             int n1 = tSizes[left].second;
             int r2 = tSizes[right].first;
             int n2 = tSizes[right].second;
            
            //coloums to copy from the left child 
             int cSI_L_1 = n1-(N-l[left].second)-1;
             int cSI_L_2 = n1-(N-l[i].second);
             
            //columns to copy from right child
             int cSI_R_1 = n1-(N-l[left].second)-1;
             int cSI_R_2 = n2-(N-l[i].second);
        }







    }    






}
