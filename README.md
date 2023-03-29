<!-- # hsseigen
Implementation of the 'divide and conquer' algorithm for Superfast HSS eigensolvers.

#To Run the code </br>
-> make sure you've installed lapack and blas. <a> http://www.netlib.org/lapack/ </a> </br>
-> In ubuntu one can also do `sudo apt-get install libblas-dev liblapack-dev liblapacke-dev` <a>https://askubuntu.com/questions/623578/installing-blas-and-lapack-packages</a></br>
-> `sudo apt install libatlas-base-dev` https://stackoverflow.com/questions/66023082/usr-bin-ld-cannot-find-ldlib-usr-bin-ld-cannot-find-lcblas-usr-bin-ld-ca <br/>
**TO RUN Programs** <br/>
-> **Make sure you've created a directory name `obj` in the project** <br/>
-> **In the terminal write the command `make` for compiling and creating all the object files** </br>
-> **after compiling the code use `./Test <filename> <matrix_size> <diagblock_size> <Band2HSS(pass 2) or Mat2HSSsym(pass 1)> <bandwidth>` to run the code**</br>
-> **by default Mat2HSSsym is used.** -->

# Hsseigen
Implementation of the 'divide and conquer' algorithm for Superfast HSS eigensolvers.

- Efficient computation of all eigenvalues of symmetric Hierarchically Semiseparable (HSS) matrices.
- We propose a shared-memory parallel algorithm and the corresponding implementation of the SuperDC eigensolver.
- Optimized the original SuperDC algorithm to reduce the storage requirement from O($N^2$) to O (N) in case of banded matrices.

### Keywords

- Symmetric eigenvalue problems
- Hierarchically semi-separable matrix
- divide-and-conquer
- Shared memory architecture

# Authors

- Anonymous

# Prerequisites

1. Install lapack and blas.

$\rightarrow$ On an ubuntu machine `sudo apt-get install libblas-dev liblapack-dev liblapacke-dev libatlas-base-dev`.

$\rightarrow$ Links to install [lapack](https://netlib.org/lapack/) and [blas](https://netlib.org/blas/) on other systems.


2. Install intel icpx compiler.

$\rightarrow$ Instructions to install icpx can be found here [intel-website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html).

3. Install [GNU make](https://www.gnu.org/software/make/) and [gcc compiler](https://gcc.gnu.org/)

# Usage 

1.  **Make sure you've created a directory name `obj` in the project.**

2. **In the terminal write the command `make` for compiling the serial version.**

3. **After compiling the code use `./Test <filename> <matrix_size> <diagblock_size> <Band2HSS(pass 2) or Mat2HSSsym(pass 1)> <bandwidth>` to run the code.**

4. **In the terminal write the command `make PARALLEL=1` for compiling the shared memory version.**

5. **After compiling the code use `./Test <filename> <matrix_size> <diagblock_size> <Band2HSS(pass 2) or Mat2HSSsym(pass 1)> <bandwidth> <no of processor>` to run the code.**


### Explanation of arguments.

| Syntax       $~~~~~~~~~~~~~~~~~~~~~~~~~$  | Description |
| :---                                      |    :----:   |
| `<filename>`                              | path to the input file containing hss matrix|
| `<matrix_size>`                           | Size of square hss matrix (can be infered from filename) |
| `<diagblock_size>`                        | Block size to split into hss structure (less than matrix size and typically 8 or 16)|
| `<Band2HSS(pass 2) or Mat2HSSsym(pass 1)>`| Whether input matrix is banded (pass 2) or not (pass 1) (can be infered from file name)|
| `<bandwidth>`                             | Half bandwidth of input matrix (can be infered from file name)|
| `<no of processor>`                       | The number of cores you want to run the shared-memory version on.|
|||
