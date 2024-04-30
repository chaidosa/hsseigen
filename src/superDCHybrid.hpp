#pragma once
// #include <superDC.h>

#ifdef HYBRD
#include "distributedRoutines.h"
#include <cmath>

SDC *HybridSuperDC(GEN *A, BinTree *btree, int *m, int mSize, int nProc, MPI_Comm process_grid)
{
    // struct timeval timeStart, timeEnd;
    // gettimeofday(&timeStart, 0);
    bt = btree;
    // cout << "Reached superDC\n";
    // Dividing Stage
    resDvd = divide2(A, bt, m, mSize);

    // cout << "Success Divide\n";
    // Conquering stage
    N = bt->GetNumNodes();

    Q0 = new EIG_MAT *[N];
    q0Sizes = new std::pair<int, int>[N];
    Lam = new double *[N];
    LamSizes = new int[N];
    l = new std::pair<int, int>[N];

    for (int k = 0; k < N; k++)
        l[k] = std::make_pair(0, 0);

    l[0] = {0, m[0] - 1};
    int it = 0;
    int lt = 0;

    for (int k = 0; k < N; k++)
    {
        std::vector<int> ch = bt->GetChildren(k + 1);
        if (ch.size() == 0)
        {
            l[k] = {lt, lt + m[it] - 1};
            lt = l[k].second + 1;
            it = it + 1;
        }
        else
        {
            l[k] = {l[ch[0] - 1].first, l[ch[1] - 1].second};
        }
    }

    for (int k = 0; k < N; k++)
        q0Sizes[k] = std::make_pair(0, 0);

    #ifdef HYBRD
    
    // Number of nodes and rank
    int nNodes, myrank;
    MPI_Comm_size(process_grid, &nNodes);
    MPI_Comm_rank(process_grid, &myrank);

    

    // int coords[2];
    // int dims[2];
    // int reorder[2];
    // MPI_Cart_get(process_grid, 2, dims, reorder, coords);
    // int my_col = coords[1]; // for row wise ordering
    // int my_row = coords[0];

    // Timer start
    double tstart, tend;

    // Divide Work btw nodes
    /**
     * Compute the leaf eigenvalues, define communicators and distribute
     */

    vector<vector<int>> level_order_nodes = bt->nodeAtLvl;
    int MaxLevel = (int)(log(nNodes)/log(2));
    vector<vector<int>> core_map(level_order_nodes.size());
    vector<vector<int>> level_map(level_order_nodes.size());
    core_map[0].push_back(0);
    level_map[0].push_back(0);
    if (level_order_nodes.size() >= 2)
    {   
        core_map[1].resize(2);
        level_map[1].resize(2);
        core_map[1][0] = 0;
        core_map[1][1] = 1;

        level_map[1][0] = 0;
        level_map[1][1] = 1;
    }

    for (int level = 2; level < level_order_nodes.size(); level++)
    {
        int num_nodes = (int)pow(2, level);
        core_map[level].resize(num_nodes);
        level_map[level].resize(num_nodes);

        int even_step = (int)pow(2, level-1) -2;
        for(int i=0; i<num_nodes/2; i++){
            if (i%2==0)
            {
                core_map[level][i] = core_map[level-1][i/2];
                level_map[level][core_map[level-1][i/2]] = i; 
            } else
            {
                core_map[level][i] = even_step+2;
                if (level > MaxLevel)
                {
                    core_map[level][i] = core_map[level][i-1];
                }
                
                level_map[level][even_step+2] = i;
                even_step+=2;
            }
        }

        int odd_step = (int)pow(2, level-1) -1;
        for (int i = num_nodes/2; i < num_nodes; i++)
        {
            if (i%2==0)
            {
                core_map[level][i] = core_map[level-1][i/2];
                level_map[level][core_map[level-1][i/2]] = i; 
            } else
            {
                core_map[level][i] = odd_step+2;
                if (level>MaxLevel){
                    core_map[level][i] = core_map[level][i-1];
                }

                level_map[level][odd_step+2] = i;
                odd_step+=2;
            }
        }
    }
    
    // Each node must solve the subtree and append to Q0_list
    tstart = MPI_Wtime();
    set<int> Q0_list;

    omp_set_num_threads(nProc);
    for (int level = level_order_nodes.size() - 1; level > MaxLevel; level--){
        
        for (int node_index = 0; node_index < level_order_nodes[level].size(); node_index++){
            if (core_map[level][node_index] == myrank)
            {   
                int node = level_order_nodes[level][node_index];
                if (bt->GetChildren(node).empty())
                {
                    int i = node - 1;
                    std::pair<double *, double *> E = computeLeafEig(make_pair(resDvd->dSizes[i].first, resDvd->dSizes[i].second), resDvd->D[i], i);

                    Lam[i] = E.first;
                    LamSizes[i] = resDvd->dSizes[i].first;
                    
                    Q0[i] = new EIG_MAT();
                    Q0[i]->Q0_leaf = E.second;
                    Q0[i]->Q0_nonleaf = NULL;
                    q0Sizes[i] = {resDvd->dSizes[i].first, resDvd->dSizes[i].second};
                    
                    Q0_list.insert(i);
                }

                else if (!bt->GetChildren(node).empty())
                {
                    vector<int> ch = bt->GetChildren(node);
                    int left = ch[0];
                    int right = ch[1];
                    int i = node - 1;
                    resDvd->Z[i] = superdcmv_desc(Q0, q0Sizes, (resDvd->Z[i]), resDvd->zSizes[i], bt, i, 1, l, fmmTrigger);

                    Lam[i] = new double[(LamSizes[left - 1]) + (LamSizes[right - 1])];
                    std::copy(Lam[left - 1], Lam[left - 1] + LamSizes[left - 1], Lam[i]);
                    std::copy(Lam[right - 1], Lam[right - 1] + LamSizes[right - 1], Lam[i] + LamSizes[left - 1]);

                    LamSizes[i] = (LamSizes[left - 1]) + (LamSizes[right - 1]);
                    delete[] Lam[left - 1];
                    delete[] Lam[right - 1];

                    LamSizes[left - 1] = 0;
                    LamSizes[right - 1] = 0;

                    int r = resDvd->zSizes[i].second;

                    Q0[i] = new EIG_MAT();
                    Q0[i]->Q0_leaf = NULL;

                    nonleaf **n_leaf = new nonleaf *[r];
                    std::pair<double *, nonleaf **> result = r_RankOneUpdate(Lam[i], LamSizes[i], resDvd->zSizes[i], resDvd->Z[i], n_leaf, r);
                    Lam[i] = result.first;
                    Q0[i]->Q0_nonleaf = result.second;
                    Q0[i]->n_non_leaf = r;
                    q0Sizes[i] = {1, r};

                    Q0_list.insert(i);
                }
            }
        }
    
    }

    for (int level = MaxLevel; level >= 0; level--)
    {   

        vector<int> work_done;
        for (int node_index = 0; node_index < level_order_nodes[level].size(); node_index++)
        {
            if (core_map[level][node_index] == myrank)
            {   
                work_done.push_back(node_index);
                int node = level_order_nodes[level][node_index];
                if (bt->GetChildren(node).empty())
                {
                    int i = node - 1;
                    std::pair<double *, double *> E = computeLeafEig(make_pair(resDvd->dSizes[i].first, resDvd->dSizes[i].second), resDvd->D[i], i);

                    Lam[i] = E.first;
                    LamSizes[i] = resDvd->dSizes[i].first;
                    
                    Q0[i] = new EIG_MAT();
                    Q0[i]->Q0_leaf = E.second;
                    Q0[i]->Q0_nonleaf = NULL;
                    q0Sizes[i] = {resDvd->dSizes[i].first, resDvd->dSizes[i].second};
                    
                    Q0_list.insert(i);
                }

                else if (!bt->GetChildren(node).empty())
                {
                    vector<int> ch = bt->GetChildren(node);
                    int left = ch[0];
                    int right = ch[1];
                    int i = node - 1;
                    resDvd->Z[i] = superdcmv_desc(Q0, q0Sizes, (resDvd->Z[i]), resDvd->zSizes[i], bt, i, 1, l, fmmTrigger);

                    Lam[i] = new double[(LamSizes[left - 1]) + (LamSizes[right - 1])];
                    std::copy(Lam[left - 1], Lam[left - 1] + LamSizes[left - 1], Lam[i]);
                    std::copy(Lam[right - 1], Lam[right - 1] + LamSizes[right - 1], Lam[i] + LamSizes[left - 1]);

                    LamSizes[i] = (LamSizes[left - 1]) + (LamSizes[right - 1]);
                    delete[] Lam[left - 1];
                    delete[] Lam[right - 1];

                    LamSizes[left - 1] = 0;
                    LamSizes[right - 1] = 0;

                    int r = resDvd->zSizes[i].second;

                    Q0[i] = new EIG_MAT();
                    Q0[i]->Q0_leaf = NULL;

                    nonleaf **n_leaf = new nonleaf *[r];
                    std::pair<double *, nonleaf **> result = r_RankOneUpdate(Lam[i], LamSizes[i], resDvd->zSizes[i], resDvd->Z[i], n_leaf, r);
                    Lam[i] = result.first;
                    Q0[i]->Q0_nonleaf = result.second;
                    Q0[i]->n_non_leaf = r;
                    q0Sizes[i] = {1, r};

                    Q0_list.insert(i);
                }
            }
        }

        int SEND_Error;
        if (myrank >= (int)pow(2, level - 1) && myrank < (int)pow(2, level) && level>0)
        { // send
            int myindx = level_map[level][myrank];
            int send_to = core_map[level][myindx - 1];

            // first send the number of results to sent
            int t = work_done.size();

            // first send contains sizes and if sending leaf or node
            SEND_Error= MPI_Send((void *)&t, 1, MPI_INT, send_to, 0, process_grid);
            if (SEND_Error!=MPI_SUCCESS)
            {
                cout << "SendError at 670" << "\n";
            }
            

            for (int i = 0; i < t; i++)
            {
                int node = level_order_nodes[level][work_done[i]];
                const int indx = node - 1;
                int LamSize = LamSizes[indx];
                int q0_size_first = q0Sizes[indx].first;
                int q0_size_second = q0Sizes[indx].second;
                bool is_leaf = bt->GetChildren(node).empty();
                double *send_buf = new double[6];

                send_buf[0] = indx;
                send_buf[1] = LamSize;
                send_buf[2] = q0_size_first;
                send_buf[3] = q0_size_second;
                send_buf[4] = is_leaf;
                send_buf[5] = Q0_list.size();

                SEND_Error= MPI_Send(send_buf, 6, MPI_DOUBLE, send_to, 0, process_grid);
                delete[] send_buf;                
                if (SEND_Error != MPI_SUCCESS)
                {
                    cout << "SendError at 696"
                         << "\n";
                }                        

                if (is_leaf)
                {

                    SEND_Error = MPI_Send(Lam[indx], LamSize, MPI_DOUBLE, send_to, 0, process_grid);
                    if (SEND_Error != MPI_SUCCESS)
                    {
                        cout << "SendError at 700"
                             << "\n";
                    }

                    SEND_Error = MPI_Send(Q0[indx]->Q0_leaf, q0_size_first * q0_size_second, MPI_DOUBLE, send_to, 0, process_grid);
                    if (SEND_Error != MPI_SUCCESS)
                    {
                        cout << "SendError at 670"
                             << "\n";
                    }
                }
                else
                {   
                    SEND_Error = MPI_Send(Lam[indx], LamSize, MPI_DOUBLE, send_to, 0, process_grid);
                    if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Lam at 700, rank: " << myrank << "\n";}

                    int* Q0_list_buff = new int[2*Q0_list.size()]; int step=0;
                    for (int x : Q0_list){
                        Q0_list_buff[step] = x;
                        step++;
                    }
                    for (int x : Q0_list){
                        Q0_list_buff[step] = Q0[x]->n_non_leaf;
                        step++;
                    }
                    SEND_Error = MPI_Send(Q0_list_buff, 2*Q0_list.size(), MPI_INT, send_to, 0, process_grid);
                    if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Q0_list at 700, rank: " << myrank << "\n";}
                    
                    

                    for (int Q0_indx : Q0_list){

                        int* Q0_indx_sizes = new int[4+Q0[Q0_indx]->n_non_leaf*(12+10+4)+2];
                        serializeQsizes(Q0[Q0_indx], Q0_indx_sizes, bt->GetChildren(Q0_indx+1).empty(), q0Sizes[Q0_indx]);
                        SEND_Error = MPI_Send(Q0_indx_sizes, 4+Q0[Q0_indx]->n_non_leaf*(12+10+4)+2, MPI_INT, send_to, 0, process_grid);
                        if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Q0_sizes at 700, rank: " << myrank << "\n";}

                        if (bt->GetChildren(Q0_indx+1).empty())
                        {
                            SEND_Error = MPI_Send(Q0[Q0_indx]->Q0_leaf, q0Sizes[Q0_indx].first*q0Sizes[Q0_indx].second, MPI_DOUBLE, send_to, 0, process_grid);
                            if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Q0_leaf at 700, rank: " << myrank << "\n";}
                        } 
                        else
                        {   
                            
                            int int_buff_sz = Q0_indx_sizes[4+Q0[Q0_indx]->n_non_leaf*(12+10+4)];
                            int double_buff_sz = Q0_indx_sizes[4+Q0[Q0_indx]->n_non_leaf*(12+10+4)+1];

                            int* T_ORG = new int[int_buff_sz];
                            double* double_buff = new double[double_buff_sz];

                            serializeQNonLeaf(Q0[Q0_indx], T_ORG, double_buff);

                            SEND_Error = MPI_Send(T_ORG, int_buff_sz, MPI_INT, send_to, 0, process_grid);
                            if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Lam at 700, rank: " << myrank << "\n";}

                            SEND_Error = MPI_Send(double_buff, double_buff_sz, MPI_DOUBLE, send_to, 0, process_grid);
                            if (SEND_Error != MPI_SUCCESS) { cout << "Cound not send Lam at 700, rank: " << myrank << "\n";}
                            // cout << "Sending sizes " << int_buff_sz << ", " << double_buff_sz << " myrank " << myrank << " sends to " << send_to << "\n";


                            // delete[] T_ORG;
                            // delete[] double_buff;
                        }

                        // delete[] Q0_indx_sizes;
                    }
                    
                    // delete[] Q0_list_buff;
                }
            }
        }
        else
        { // recv

            if (myrank < (int)pow(2, level) && level>0)
            {

                int myindx = level_map[level][myrank];
                int recv_from = core_map[level][myindx + 1];

                // number of results to recv
                int t = 0;
                MPI_Status status;
                MPI_Recv((void *)&t, 1, MPI_INT, recv_from, 0, process_grid, &status);

                for (int i = 0; i < t; i++)
                {
                    double *recv_sizes_first = new double[6];
                    MPI_Recv((void *)recv_sizes_first, 6, MPI_DOUBLE, recv_from, 0, process_grid, &status);

                    int indx = (int)recv_sizes_first[0];
                    int lamsz = (int)recv_sizes_first[1];
                    int qsz_first = (int)recv_sizes_first[2];
                    int qsz_second = (int)recv_sizes_first[3];
                    int isLeaf = (int)recv_sizes_first[4];
                    int Q0_list_size = (int)recv_sizes_first[5];


                    q0Sizes[indx] = make_pair(qsz_first, qsz_second);
                    LamSizes[indx] = lamsz;

                    // cout << "Myrank, recvd Qi: " << myrank << "," << indx << "\n";
                    if (isLeaf)
                    {
                        Lam[indx] = new double[lamsz];
                        Q0[indx] = new EIG_MAT();
                        Q0[indx]->Q0_nonleaf = NULL;
                        Q0[indx]->Q0_leaf = new double[qsz_first * qsz_first];
                        MPI_Recv(Lam[indx], lamsz, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                        MPI_Recv(Q0[indx]->Q0_leaf, qsz_first * qsz_first, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                        Q0_list.insert(indx);
                    }
                    else
                    {

                        Lam[indx] = new double[lamsz];
                        MPI_Recv(Lam[indx], lamsz, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                        

                        int* Q0_list_buff = new int[2*Q0_list_size];
                        MPI_Recv(Q0_list_buff, 2*Q0_list_size, MPI_INT, recv_from, 0, MPI_COMM_WORLD, &status);

                        

                        for (int _indx = 0; _indx < Q0_list_size; _indx++)
                        {   
                            int step = Q0_list_size;
                            int Q0_n_non_leaf = Q0_list_buff[step+_indx];
                            int Q0_indx = Q0_list_buff[_indx];
                            Q0[Q0_indx] = new EIG_MAT();
                            Q0[Q0_indx]->n_non_leaf = Q0_n_non_leaf;
                            int* Q0_sizes_buff = new int[4+Q0_n_non_leaf*(12+10+4)+2]; // error here
                            MPI_Recv(Q0_sizes_buff, 4+Q0_n_non_leaf*(12+10+4)+2, MPI_INT, recv_from, 0, process_grid, &status);
                            bool isleaf = Q0_sizes_buff[2];


                            

                            pair<int,int> buff_sizes = deserializeQsizes(Q0[Q0_indx], Q0_sizes_buff);
                            q0Sizes[Q0_indx] = make_pair(Q0_sizes_buff[0], Q0_sizes_buff[1]);
                            if (isleaf)
                            {
                                Q0[Q0_indx]->Q0_leaf = new double[q0Sizes[Q0_indx].first*q0Sizes[Q0_indx].second];
                                MPI_Recv(Q0[Q0_indx]->Q0_leaf, q0Sizes[Q0_indx].first*q0Sizes[Q0_indx].second, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                            }
                            else
                            {
                                int int_buff_sz = buff_sizes.first;
                                int double_buff_sz = buff_sizes.second;

                                int* T_ORG = new int[int_buff_sz];
                                double* double_buff = new double[double_buff_sz];

                                MPI_Recv(T_ORG, int_buff_sz, MPI_INT, recv_from, 0, process_grid, &status);

                                MPI_Recv(double_buff, double_buff_sz, MPI_DOUBLE, recv_from, 0, process_grid, &status);
                                // if (level==1&&myrank==0)
                                // cout << "Reached Here 825 " << myrank << " " << Q0_indx << " " << Q0_sizes_buff[2] << " " << Q0_sizes_buff[3] << " " << qsz_second << "\n";
                                deserialzeQNonLeaf(Q0[Q0_indx], T_ORG, double_buff);
                                
                            }

                            Q0_list.insert(Q0_indx);
                            delete[] Q0_sizes_buff;
                        }
                        
                        delete[] Q0_list_buff;
                    }
                }
            }
        }

        // MPI_Barrier(process_grid);
        // cout << "Completed level " << level << " " << myrank << "\n";
        // if(myrank==0){
        //     cout << "\n";
        // }
        // sleep(1);
        MPI_Barrier(process_grid);
    }

    

    tend = MPI_Wtime();
    double etime = tend - tstart;
    double max_etime = 0.0;
    MPI_Reduce(&etime, &max_etime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0)
    {
        cout << "Hybrid Superdc took " << max_etime << " seconds"
             << "\n";
    }
    MPI_Barrier(process_grid);
    
    #endif

    #ifdef HYBRD
    if (myrank==0){
    #endif
    std::ofstream txtOut;
    txtOut.open("output.txt", std::ofstream::out);
    // txtOut << setprecision(10) << elapsed / (double)1000000 << " seconds" << endl;
    vector<double> tempeig;
    for (int k = 0; k < LamSizes[N - 1]; k++)
        tempeig.push_back(Lam[N - 1][k]);

    std::sort(tempeig.begin(), tempeig.end());
    int count = 0;
    for (int k = 0; k < LamSizes[N - 1]; k++)
    {
        count++;
        txtOut << setprecision(20) << tempeig[k] << endl;
    }
    cout << count;
    SDC *resSDC = new SDC();

    resSDC->Q = Q0;
    resSDC->L = Lam[N - 1];
    resSDC->lSize = LamSizes[N - 1];
    resSDC->qSizes = q0Sizes;
    return resSDC;
    
    #ifdef HYBRD
    } 
    else
    {
        return new SDC();
    }
    #endif

    return nullptr;
}
#endif
