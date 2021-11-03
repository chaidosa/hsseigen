#include<string.h>
#include<algorithm>
#include<bitset>
#include<assert.h>
#include<stdio.h>
#include"QR.h"
extern "C"
{
    #include<lapacke.h>
    #include<lapack.h>
}

void qpr( double* const _Q, double* const _R, int* const _P, double* const _A, const size_t _m, const size_t _n) {
    // Maximal rank is used by Lapacke
    const size_t rank = std::min(_m, _n); 

    // Tmp Array for Lapacke
    double* tau = new double[rank];

    // Calculate QR factorisations
    lapack_int info = LAPACKE_dgeqpf(LAPACK_ROW_MAJOR, (int) _m, (int) _n, _A, (int) _n, _P, tau);
    if(info != 0)
    {
	printf("ERROR %d parameter had an illegal value\n",info);
	assert(0);
    }

    // Copy the upper triangular Matrix R (rank x _n) into position
    for(size_t row =0; row < rank; ++row) {
        memset(_R+row*_n, 0, row*sizeof(double)); // Set starting zeros
        memcpy(_R+row*_n+row, _A+row*_n+row, (_n-row)*sizeof(double)); // Copy upper triangular part from Lapack result.
    }

    // Create orthogonal matrix Q (in tmpA)
    info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, (int) _m, (int) rank, (int) rank, _A, (int) _n, tau);
    if(info != 0)
    {
	printf("%d parameter had an illegal value\n",info);
	assert(0);
    }

    //Copy Q (_m x rank) into position
    if(_m == _n) {
        memcpy(_Q, _A, sizeof(double)*(_m*_n));
    } else {
        for(size_t row =0; row < _m; ++row) {
            memcpy(_Q+row*rank, _A+row*_n, sizeof(double)*(rank));
        }
    }
	
    delete [] tau;
}


//source: http://www.geeksforgeeks.org/inplace-m-x-n-size-matrix-transpose/
void GetTransposeInPlace(double *A, int r, int c)
{
    using namespace std;
    int N = r*c;
    int size = N - 1;
    double t; // holds element to be replaced, eventually becomes next element to move
    int next; // location of 't' to be moved
    int cycleBegin; // holds start of cycle
    int i; // iterator
    string b(N,0); //bitset.
 
    b[0] = b[N-1] = 1;
    i = 1; // Note that A[0] and A[N-1] won't move
    while (i < size)
    {
        cycleBegin = i;
        t = A[i];
        do
        {
            // i_new = (i*r)%(N-1)
            next = (i*r)%size;
            swap(A[next], t);
            b[i] = 1;
            i = next;
        }
        while (i != cycleBegin);
 
        for (i = 1; i < size && b[i]; i++);
    }
}
