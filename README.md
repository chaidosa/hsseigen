# Hsseigen

Implementation of the 'divide and conquer' algorithm for Superfast HSS eigensolvers.

## Overview

- Efficient computation of all eigenvalues of symmetric Hierarchically Semiseparable (HSS) matrices.
- A shared-memory parallel algorithm and the corresponding implementation of the SuperDC eigensolver.
- Optimized the original SuperDC algorithm to reduce storage requirements from O($N^2$) to O(N) for banded matrices.

📄 **Paper:** [ACM Digital Library](https://dl.acm.org/doi/10.1145/3673038.3673119)

### Keywords
- Symmetric eigenvalue problems
- Hierarchically semi-separable matrix (HSS matrix)
- Divide-and-conquer
- Shared memory architecture

## Authors
Abhishek Josyula, Pritesh Verma, Amar Gaonkar, Amlan Barua, Nikhil Hegde

---

## Prerequisites

### 1. Install LAPACK, BLAS, or OpenBLAS

- **Ubuntu:**
  ```sh
  sudo apt-get install libblas-dev liblapack-dev liblapacke-dev libatlas-base-dev libopenblas-dev
  ```
- **Other OS:**
  - [LAPACK](https://netlib.org/lapack/)
  - [BLAS](https://netlib.org/blas/)
  - [OpenBLAS](https://www.openblas.net/)

### 2. Install Intel ICPX Compiler
- Instructions available on the [Intel website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html).

### 3. Install GNU Make and GCC Compiler
- [GNU Make](https://www.gnu.org/software/make/)
- [GCC Compiler](https://gcc.gnu.org/)

---

## CodeSpaces Support

1. Fork the repository and open it in **GitHub CodeSpaces** to work with all prerequisites installed.
2. Click on the **`<> Code`** button.
3. Choose **"Create CodeSpace on master"** to launch the development environment.

---

## Usage

### Compilation

- **Serial Version:**
  ```sh
  make
  ```
- **Parallel Version (Shared Memory):**
  ```sh
  make PARALLEL=1
  ```
- **Binary Output:**
  - The compiled binary `Test` will be created in the root directory.

### Execution

```sh
./Test <filename> <matrix_size> <diagblock_size> <Band2HSS(pass 2) or Mat2HSSsym(pass 1)> <bandwidth>
```

**For parallel execution:**
```sh
./Test <filename> <matrix_size> <diagblock_size> <Band2HSS(pass 2) or Mat2HSSsym(pass 1)> <bandwidth> <no_of_processors>
```

### Argument Explanation

| Argument | Description |
|----------|-------------|
| `<filename>` | Path to the input file containing the HSS matrix. |
| `<matrix_size>` | Size of the square HSS matrix (can be inferred from the filename). |
| `<diagblock_size>` | Block size to split into HSS structure (should be less than matrix size, typically 8 or 16). |
| `<Band2HSS(pass 2) or Mat2HSSsym(pass 1)>` | Specifies whether the input matrix is banded (pass 2) or not (pass 1) (can be inferred from the filename). |
| `<bandwidth>` | Half bandwidth of the input matrix (can be inferred from the filename). |
| `<no_of_processors>` | Number of cores to use in the shared-memory version. |

---

## Input Data

📂 **Download HSS structured test matrices:** [Google Drive Link](https://drive.google.com/drive/folders/1Qs-U8bQf_apAt8LKpTy5IX-W5Y0RBPGk)

---

## Miscellaneous

- The project has been tested using OpenBLAS, LAPACK, and Intel toolchains.

