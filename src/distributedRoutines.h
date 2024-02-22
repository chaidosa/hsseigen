#ifndef MPI_HELPER_FUNC
#define MPI_HELPER_FUNC

#include<bits/stdc++.h>
#include<string.h>
#include "BinTree.h"
#include "superDC.h"
#include "divide2.h"
#include "superdcmv_desc.h"
#include "superdcmv_cauchy.h"
#include "eigenmatrix.h"
#include "secular.h"
#include "Generators.h"
#include <sys/time.h>

#ifdef DIST
// extern "C"{
#include <mpi.h>
// }
#endif

void printHssTree(vector<vector<int>> level_order_nodes, vector<vector<int>> core_map);

void serializeQNonLeaf(EIG_MAT* serialze, int* T_ORG, double* double_buff);

void serializeQsizes(const EIG_MAT* serialize, int* buff, bool isleaf, pair<int,int> q0size);

pair<int,int> deserializeQsizes(EIG_MAT* deserialize, int* buff);

void deserialzeQNonLeaf(EIG_MAT* deserialze, int* T_ORG, double* double_buff);
#endif