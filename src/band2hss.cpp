#include<string.h>
#include<math.h>
#include<algorithm>
#include<assert.h>
#include "BinTree.h"
#include "band2hss.h"
#include "Generators.h"

#if defined(DIST) || defined(HYBRD) || defined(PARALLEL_TASK) || defined(OPEN_CILK)
GEN *band2hss(double *AA, int aSize, BinTree* bt, int* m, int mSize, int w){
    /**
     * @brief Computes HSS form of a banded matrix A
     * @param    AA: banded form
     *           w:  half bandwidth, total bandwidth = 2*w+1
     *           bt: postordered binary tree
     *           m:  partition, e.g., m = [100 100 100 100] when n = 400
     * 
     * @return  D, U, R, B Generators:HSS of the banded form of Class B2HSS
     */
    
    GEN * result = new GEN();
    double *A = AA;
    int n = bt->GetNumNodes();

    int aRowWidth         = sqrt(aSize);
    int aColWidth         = aRowWidth;

    int numElemsInaRow = aRowWidth;
    int maxElemsInaColumn = w+w+1;
     //To store the diagonal blocks
    double** D = new double*[n];

    //create a list to hold the matrix dimensions of the diagonal matrices 
	std::pair<int, int>* dSizes = new std::pair<int, int>[n];
    for(int i = 0; i < n; i++)
		dSizes[i]=std::make_pair(0,0);

     //To store the Generator U
    double** U = new double*[n];

    //create a list to hold the matrix dimensions of the generator U
	std::pair<int, int>* uSizes = new std::pair<int, int>[n];
    for(int i = 0; i < n; i++)
		uSizes[i]=std::make_pair(0,0);

    //to store generator R
    double** R = new double*[n];

    //create a list to hold the matrix dimensions of the generator R
	std::pair<int, int>* rSizes = new std::pair<int, int>[n];
    for(int i = 0; i < n; i++)
		rSizes[i]=std::make_pair(0,0);

    //to store generator B
    double** B = new double*[n];

    //create a list to hold the matrix dimensions of the generator B
	std::pair<int, int>* bSizes = new std::pair<int, int>[n];

    for(unsigned int i = 0; i < n; i++)
		bSizes[i]=std::make_pair(0,0);
    if(n == 1){
        // Return type to be defined here
        return NULL;
    }

    vector<double> cl;
    for(unsigned int i = 1; i<=n; i++){
        std::vector<int> ch = bt->GetChildren(i);
        
        // if the node is child node
        if(ch.size() == 0)
            cl.push_back(i);
    }

    vector<double> l;
    l.push_back(0);
    for(unsigned int i = 0; i < mSize; ++i){
        double temp = l[i] + m[i];
        l.push_back(temp);
    }

    //Zeros
    double *O = new double[w*w];
    memset(O, 0, sizeof(double)*(w*w));
    std::pair<int, int> OSize = {w, w};
    
    //identity matrix
    double *I = new double[w];
    memset(I, 0, sizeof(double)*w);

    int k = -1;

    for(unsigned int i = 1; i < n; i++)
    {
        std::vector<int> ch = bt->GetChildren(i);
        //if the node in child node
        if(ch.size() == 0)
        {
            k += 1;
            int dRow = l[k+1]-l[k];
            int dCol = dRow;

            dSizes[i-1] = {dRow, dCol};
            D[i-1] = new double[dRow*dCol];
            double *tempD = D[i-1];
            memset(tempD, 0, sizeof(double)*dRow*dCol);
            // copying the diagonal elements from A
            //int D_filled = 0;
            for(unsigned int j = (int)l[k], D_filled = 0; j < (int)(l[k+1]); j++, D_filled++)
            {
                
                for(int col = 0; col < dCol; col++){
                    if((j - (col+(int)l[k]) + w) < 0 || (j - (col+(int)l[k]) + w) >= maxElemsInaColumn){
                        tempD[col + D_filled*dCol] = 0;
                       // continue;
                    }
                    else{
                        int loc = col + (int)l[k] + (j- (col + (int)l[k]) + w)*numElemsInaRow;
                        tempD[col + D_filled*dCol] = A[loc];
                    }
                }
                
                //memcpy(tempD + D_filled, A+((int)l[k] +j*aRowWidth), sizeof(double)*dCol);
                //D_filled++;
            }

            int uRow = m[k];
            int uCol = 2*w;
            uSizes[i-1] = {uRow, uCol};

            U[i-1] = new double[uRow*uCol];
            memset(U[i-1], 0, sizeof(double)*(uRow*uCol));

            for(unsigned int j = 0; j < w; j++)
                U[i-1][j + j*uCol] = 1;

            int startROw = uRow - w; //for our case indexing starts from 0;
            int startCol = w;
        
            for(unsigned int row = startROw,col = startCol; row < uRow, col < 2*w; row++, col++)
            {
                    U[i-1][col + row*uCol] = 1;
            }
        }
        int parentID = bt->tr[i-1];
        int rc = 0;
        int ii;
        ch.clear();
        ch = bt->GetChildren(parentID);
        if(ch[0] == i){
            rc = i;

            while(bt->GetChildren(rc).size() != 0)
            {
                rc = bt->GetChildren(rc)[1];
            } 

            for(unsigned int iter = 0; iter < cl.size(); iter++)
            {
                if((int)cl[iter] == rc)
                {
                    ii = iter;
                    break;
                }
            }

            int bRow = 2*w;
            int bCol = bRow;
            bSizes[i-1] = {bRow, bCol};
            B[i-1] = new double[bRow*bCol];
            memset(B[i-1], 0, sizeof(double)*(bRow*bCol));
        
            int aRowStart = l[ii+1]-w;
            int aRowEnd   = l[ii+1];
            int aColStart = l[ii+1];
            int aColEnd   = l[ii+1]+w-1; 
            int iter = aColStart;
            for(unsigned int row = aRowStart, t = w; row < aRowEnd; row++,t++){
                for(int col = 0, iter = aColStart; col < w; col++, iter++){
                    if((row - iter + w) < 0 || (row - iter + w) >= maxElemsInaColumn){
                        B[i-1][col + t*bCol] = 0;
                    }
                    else{
                        B[i-1][col + t*bCol] = A[iter + (row - iter + w)*numElemsInaRow];
                    }
                }               
            }
                //memcpy(B[i-1]+(t*bCol), A+aColStart+(row*aRowWidth), sizeof(double)*(w));
        
            int rCol = 2*w;
            int rRow = rCol;
            R[i-1] = new double[(rRow*rCol)];
            rSizes[i-1] = {rRow, rCol};

            memset(R[i-1], 0, sizeof(double)*(rCol*rRow));

            for(unsigned int row = 0; row < w; row++)
            {
                R[i-1][row + row*rCol] = 1;
            }

        }
        else
        {
            rc = (bt->GetChildren(bt->tr[i-1]))[0];
            while(bt->GetChildren(rc).size() != 0)
            {
                rc = bt->GetChildren(rc)[1];
            }

            for(unsigned int iter = 0; iter < cl.size(); iter++)
            {
                if((int)cl[iter] == rc)
                {
                    ii = iter;
                    break;
                }
            }

            int bRow = 2*w;
            int bCol = bRow;
            bSizes[i-1] = {bRow, bCol};
            B[i-1] = new double[bRow*bCol];
            memset(B[i-1], 0, sizeof(double)*(bRow*bCol));

            int aRowStart = l[ii+1];
            int aRowEnd   = l[ii+1]+w;
            int aColStart = l[ii+1]-w;
            int aColEnd   = l[ii+1]; 
            int iterI = aRowStart;
            int iterJ = aColStart;
            for(unsigned int row = 0, iterI = aRowStart; row < w; row++, iterI++)
            {
               
                for(int col = w, iterJ = aColStart; col < 2*w; col++, iterJ++){
                    if((iterI - iterJ + w) < 0 || (iterI - iterJ + w) >= maxElemsInaColumn){
                        B[i-1][col + row*bCol] = 0;
                    }
                    else{
                        B[i-1][col + row*bCol] = A[iterJ + (iterI - iterJ + w)*numElemsInaRow];
                    }
            
               }
               // mempcpy(B[i-1]+(w)+(row*bCol), A+aColStart+(aRowStart*aRowWidth), sizeof(double)*w);
               // aRowStart++;
            }
            int rRow = 2*w;
            int rCol = rRow;
            R[i-1]   = new double[rRow*rCol];
            rSizes[i-1] = {rRow, rCol};
            memset(R[i-1], 0, sizeof(double)*(rRow*rCol));
            int col = w;
            for(int row = w; row < 2*w; row++){
                R[i-1][col + row*rCol] = 1;
                ++col;
            }
        }
    }
    result->D = D;
    result->U = U;
    result->R = R;
    result->B = B;
    result->bSizes = bSizes;
    result->rSizes = rSizes;
    result->uSizes = uSizes;
    result->dSizes = dSizes;

    return result;
}


#else
GEN *band2hss(double *AA, int aSize, BinTree* bt, int* m, int mSize, int w){
    /**
     * @brief Computes HSS form of a banded matrix A
     * @param    AA: banded form
     *           w:  half bandwidth, total bandwidth = 2*w+1
     *           bt: postordered binary tree
     *           m:  partition, e.g., m = [100 100 100 100] when n = 400
     * 
     * @return  D, U, R, B Generators:HSS of the banded form of Class B2HSS
     */
    
    GEN * result = new GEN();
    double *A = AA;
    int n = bt->GetNumNodes();

    int aRowWidth         = sqrt(aSize);
    int aColWidth         = aRowWidth;

     //To store the diagonal blocks
    double** D = new double*[n];
    for (int i = 0; i < n; i++) D[i] = NULL;
    
    

    //create a list to hold the matrix dimensions of the diagonal matrices 
	std::pair<int, int>* dSizes = new std::pair<int, int>[n];
    for(int i = 0; i < n; i++)
		dSizes[i]=std::make_pair(0,0); 

     //To store the Generator U
    double** U = new double*[n];
    for (int i = 0; i < n; i++) U[i] = NULL;

    //create a list to hold the matrix dimensions of the generator U
	std::pair<int, int>* uSizes = new std::pair<int, int>[n];
    for(int i = 0; i < n; i++)
		uSizes[i]=std::make_pair(0,0);   

    //to store generator R
    double** R = new double*[n];
    for (int i = 0; i < n; i++) R[i] = NULL;


    //create a list to hold the matrix dimensions of the generator R
	std::pair<int, int>* rSizes = new std::pair<int, int>[n];
    for(int i = 0; i < n; i++)
		rSizes[i]=std::make_pair(0,0);

    //to store generator B
    double** B = new double*[n];
    for (int i = 0; i < n; i++) B[i] = NULL;


    //create a list to hold the matrix dimensions of the generator B
	std::pair<int, int>* bSizes = new std::pair<int, int>[n];

    for(unsigned int i = 0; i < n; i++)
		bSizes[i]=std::make_pair(0,0);                      

    if(n == 1){
        // Return type to be defined here
        return NULL;
    }

    vector<double> cl;
    for(unsigned int i = 1; i<=n; i++){
        std::vector<int> ch = bt->GetChildren(i);
        
        // if the node is child node
        if(ch.size() == 0)
            cl.push_back(i);      

    }

    vector<double> l;
    l.push_back(0);
    for(unsigned int i = 0; i < mSize; ++i){
        double temp = l[i] + m[i];
        l.push_back(temp);
    }

    //Zeros
    double *O = new double[w*w];
    memset(O, 0, sizeof(double)*(w*w));
    std::pair<int, int> OSize = {w, w};
    
    //identity matrix
    double *I = new double[w];
    memset(I, 0, sizeof(double)*w);

    int k = -1;

    for(unsigned int i = 1; i < n; i++)
    {
        std::vector<int> ch = bt->GetChildren(i);
        //if the node in child node
        if(ch.size() == 0)
        {
            k += 1;
            int dRow = l[k+1]-l[k];
            int dCol = dRow;

            dSizes[i-1] = {dRow, dCol};
            D[i-1] = new double[dRow*dCol];
            double *tempD = D[i-1];
            // copying the diagonal elements from A
            int D_filled = 0;
            for(unsigned int j = l[k]; j < (l[k+1]); j++)
            {
                memcpy(tempD + D_filled, A+((int)l[k] +j*aRowWidth), sizeof(double)*dCol);
                D_filled += dCol;
            }

            int uRow = m[k];
            int uCol = 2*w;
            uSizes[i-1] = {uRow, uCol};

            U[i-1] = new double[uRow*uCol];
            memset(U[i-1], 0, sizeof(double)*(uRow*uCol));

            for(unsigned int j = 0; j < w; j++)
                U[i-1][j + j*uCol] = 1;

            int startROw = uRow - w; //for our case indexing starts from 0;
            int startCol = w;
        
            for(unsigned int row = startROw,col = startCol; row < uRow, col < 2*w; row++, col++)
            {                                    
                    U[i-1][col + row*uCol] = 1;                    
            }          
        }        

        int parentID = bt->tr[i-1];
        int rc = 0;
        int ii;
        ch.clear();
        ch = bt->GetChildren(parentID);
        if(ch[0] == i){
            rc = i;

            while(bt->GetChildren(rc).size() != 0)
            {
                rc = bt->GetChildren(rc)[1];
            } 

            for(unsigned int iter = 0; iter < cl.size(); iter++)
            {
                if((int)cl[iter] == rc)
                {
                    ii = iter;
                    break;
                }
            }

            int bRow = 2*w;
            int bCol = bRow;
            bSizes[i-1] = {bRow, bCol};
            B[i-1] = new double[bRow*bCol];
            memset(B[i-1], 0, sizeof(double)*(bRow*bCol));
        
            int aRowStart = l[ii+1]-w;            
            int aRowEnd   = l[ii+1];
            int aColStart = l[ii+1];
            int aColEnd   = l[ii+1]+w-1; 
        
            for(unsigned int row = aRowStart, t = w; row < aRowEnd; row++,t++)
                memcpy(B[i-1]+(t*bCol), A+aColStart+(row*aRowWidth), sizeof(double)*(w));
        
            int rCol = 2*w;
            int rRow = rCol;
            R[i-1] = new double[(rRow*rCol)];
            rSizes[i-1] = {rRow, rCol};

            memset(R[i-1], 0, sizeof(double)*(rCol*rRow));

            for(unsigned int row = 0; row < w; row++)
            {
                R[i-1][row + row*rCol] = 1;
            }

        }
        else
        {
            rc = (bt->GetChildren(bt->tr[i-1]))[0];
            while(bt->GetChildren(rc).size() != 0)
            {
                rc = bt->GetChildren(rc)[1];
            }

            for(unsigned int iter = 0; iter < cl.size(); iter++)
            {
                if((int)cl[iter] == rc)
                {
                    ii = iter;
                    break;
                }
            }

            int bRow = 2*w;
            int bCol = bRow;
            bSizes[i-1] = {bRow, bCol};
            B[i-1] = new double[bRow*bCol];
            memset(B[i-1], 0, sizeof(double)*(bRow*bCol));

            int aRowStart = l[ii+1];
            int aRowEnd   = l[ii+1]+w;
            int aColStart = l[ii+1]-w;            
            int aColEnd   = l[ii+1]; 

            for(unsigned int row = 0; row < w; row++)
            {
                mempcpy(B[i-1]+(w)+(row*bCol), A+aColStart+(aRowStart*aRowWidth), sizeof(double)*w);
                aRowStart++;
            }
            int rRow = 2*w;
            int rCol = rRow;
            R[i-1]   = new double[rRow*rCol];
            rSizes[i-1] = {rRow, rCol};
            memset(R[i-1], 0, sizeof(double)*(rRow*rCol));
            int col = w;
            for(int row = w; row < 2*w; row++){
                R[i-1][col + row*rCol] = 1;
                ++col;
            }              
        }
    }
    result->D = D;
    result->U = U;
    result->R = R;
    result->B = B;
    result->bSizes = bSizes;
    result->rSizes = rSizes;
    result->uSizes = uSizes;
    result->dSizes = dSizes;

    return result;    
}
#endif