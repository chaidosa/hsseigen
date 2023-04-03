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
- A shared-memory parallel algorithm and the corresponding implementation of the SuperDC eigensolver.
- Optimized the original SuperDC algorithm to reduce the storage requirement from O($N^2$) to O (N) in case of banded matrices.

### Keywords

- Symmetric eigenvalue problems.
- Hierarchically semi-separable matrix or HSS matrix.
- Divide-and-conquer.
- Shared memory architecture.

# Authors

- Anonymous

# Prerequisites

1. Install lapack, blas OR openblas.

$\rightarrow$ On an ubuntu machine `sudo apt-get install libblas-dev liblapack-dev liblapacke-dev libatlas-base-dev libopenblas-dev`.

$\rightarrow$ Links to install [lapack](https://netlib.org/lapack/), [blas](https://netlib.org/blas/) and [openblas](https://www.openblas.net/) on other systems.

2. Install intel icpx compiler.

$\rightarrow$ Instructions to install icpx can be found here [intel-website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html).

3. Install [GNU make](https://www.gnu.org/software/make/) and [gcc compiler](https://gcc.gnu.org/).

# CodeSpaces Support

1. To begin working on the project with all the pre-requsites satisfied, fork the repository and begin working in Code Spaces.

2. After forking the repository, click on the **`<> Code`** button.

3. You should see option to create code space on master.

4. Click on the **`Create code space on master`** button and you will be redirected to a working environment.

# Usage 

## Compilation

1. **In the terminal write the command `make` for compiling the serial version.**

2. **In the terminal write the command `make PARALLEL=1` for compiling the shared memory version.**

3. **A binary called `Test` is created in the root directory of the project.**


## How to run
1. **After compiling the code use `./Test <filename> <matrix_size> <diagblock_size> <Band2HSS(pass 2) or Mat2HSSsym(pass 1)> <bandwidth>` to run the code.**

2. **After compiling the code use `./Test <filename> <matrix_size> <diagblock_size> <Band2HSS(pass 2) or Mat2HSSsym(pass 1)> <bandwidth> <no of processor>` to run the code.**


### Explanation of arguments.

| Syntax       $~~~~~~~~~~~~~~~~~~~~~~~~~$  | Description |
| :---                                      |    :----:   |
| `<filename>`                              | path to the input file containing hss matrix|
| `<matrix_size>`                           | Size of square hss matrix (can be infered from filename) |
| `<diagblock_size>`                        | Block size to split into hss structure (less than matrix size and typically 8 or 16)|
| `<Band2HSS(pass 2) or Mat2HSSsym(pass 1)>`| Whether input matrix is banded (pass 2) or not (pass 1) (can be infered from file name)|
| `<bandwidth>`                             | Half bandwidth of input matrix (can be infered from file name)|
| `<no of processor>`                       | The number of cores you want to run the shared-memory version on.|



# Drive link to input data:

1. You can download HSS structured test matrices [here](https://drive.google.com/drive/folders/1Qs-U8bQf_apAt8LKpTy5IX-W5Y0RBPGk)

# Miscellaneous

1. We have tested the project using openblas, lapack, and intel toolchains
