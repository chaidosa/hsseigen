#include "distributedRoutines.h"
#include <string.h>
#ifdef DIST
void printHssTree(vector<vector<int>> level_order_nodes, vector<vector<int>> core_map)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0)
    {
        for (int i = 0; i < level_order_nodes.size(); i++)
        {
            cout << "At level: " << i << " nodes are: ";
            for (int j = 0; j < level_order_nodes[i].size(); j++)
            {
                cout << level_order_nodes[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n\n";
        for (int i = 0; i < core_map.size(); i++)
        {
            cout << "At level: " << i << " nodes are computed by ranks: ";
            for (int j = 0; j < core_map[i].size(); j++)
            {
                cout << core_map[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n\n";
    }
}
#endif

void serializeQsizes(const EIG_MAT* serialize, int* buff, bool isleaf, pair<int,int> q0size){
    buff[0] = q0size.first;
    buff[1] = q0size.second;
    buff[2] = isleaf;

    if (isleaf)
    {
        buff[3] = serialize->n_non_leaf;
    }
    else
    {
        buff[3] = serialize->n_non_leaf;
    }

    int int_buff_sz = 0; 
    int double_buff_sz = 0;

    if (!isleaf)
    {
        for (int i = 0; i < serialize->n_non_leaf; i++)
        {
            int buff_indx = (12 + 10 + 4) * i + 4;
            int step = 0;
            for (int k = 0; k < 6; k++)
            {
                buff[buff_indx + step] = serialize->Q0_nonleaf[i]->qcSizes[k].first;
                step += 1;
                buff[buff_indx + step] = serialize->Q0_nonleaf[i]->qcSizes[k].second;
                step += 1;

                if (k==5)
                {
                    int_buff_sz+=serialize->Q0_nonleaf[i]->qcSizes[k].first * serialize->Q0_nonleaf[i]->qcSizes[k].second;
                } else {
                    double_buff_sz+=serialize->Q0_nonleaf[i]->qcSizes[k].first * serialize->Q0_nonleaf[i]->qcSizes[k].second;
                }
                
            }

            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->JSize.first;
            step += 1;
            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->JSize.second;
            step += 1;
            double_buff_sz+=serialize->Q0_nonleaf[i]->JSize.first*serialize->Q0_nonleaf[i]->JSize.second;

            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->GSize.first;
            step += 1;
            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->GSize.second;
            step += 1;
            double_buff_sz+=serialize->Q0_nonleaf[i]->GSize.first*serialize->Q0_nonleaf[i]->GSize.second;


            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->ISize.first;
            step += 1;
            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->ISize.second;
            step += 1;
            double_buff_sz+=serialize->Q0_nonleaf[i]->ISize.first*serialize->Q0_nonleaf[i]->ISize.second;


            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->v2cSize.first;
            step += 1;
            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->v2cSize.second;
            step += 1;
            double_buff_sz+=serialize->Q0_nonleaf[i]->v2cSize.first*serialize->Q0_nonleaf[i]->v2cSize.second;


            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->TSize.first;
            step += 1;
            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->TSize.second;
            step += 1;
            int_buff_sz+=serialize->Q0_nonleaf[i]->TSize.first*serialize->Q0_nonleaf[i]->TSize.second;

            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->n;
            step += 1;
            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->n1;
            step += 1;
            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->n2;
            step += 1;
            buff[buff_indx + step] = serialize->Q0_nonleaf[i]->n3;
            step += 1;
        }

        buff[4 + serialize->n_non_leaf * (12 + 10 + 4)] = int_buff_sz;
        buff[4 + serialize->n_non_leaf * (12 + 10 + 4) + 1] = double_buff_sz;
    }
}

pair<int,int> deserializeQsizes(EIG_MAT* deserialize, int* buff){

    bool isleaf = buff[2];
    // deserialize->n_non_leaf = buff[3];
    
    if (!isleaf)
    {   
        
        deserialize->Q0_nonleaf = new nonleaf*[deserialize->n_non_leaf];
        
        for (int i = 0; i < deserialize->n_non_leaf; i++)
        {
            int buff_indx = (12+10+4)*i+4;
            int step = 0;

            deserialize->Q0_nonleaf[i] = new nonleaf();

            for (int j = 0; j < 6; j++)
            {
                deserialize->Q0_nonleaf[i]->qcSizes[j] = make_pair(buff[buff_indx+step], buff[buff_indx+step+1]);
                step+=2;
            }

            deserialize->Q0_nonleaf[i]->JSize = make_pair(buff[buff_indx + step], buff[buff_indx + step + 1]);
            step += 2;

            deserialize->Q0_nonleaf[i]->GSize = make_pair(buff[buff_indx + step], buff[buff_indx + step + 1]);
            step += 2;

            deserialize->Q0_nonleaf[i]->ISize = make_pair(buff[buff_indx + step], buff[buff_indx + step + 1]);
            step += 2;

            deserialize->Q0_nonleaf[i]->v2cSize = make_pair(buff[buff_indx + step], buff[buff_indx + step + 1]);
            step += 2;

            deserialize->Q0_nonleaf[i]->TSize = make_pair(buff[buff_indx + step], buff[buff_indx + step + 1]);
            step += 2;

            deserialize->Q0_nonleaf[i]->n = buff[buff_indx+step];
            step+=1;

            deserialize->Q0_nonleaf[i]->n1 = buff[buff_indx+step];
            step+=1;

            deserialize->Q0_nonleaf[i]->n2 = buff[buff_indx+step];
            step+=1;

            deserialize->Q0_nonleaf[i]->n3 = buff[buff_indx+step];
            step+=1;
            
        }
        
    }
    
    return make_pair(buff[4+deserialize->n_non_leaf*(12+10+4)], buff[4+deserialize->n_non_leaf*(12+10+4)+1]);
}


void serializeQNonLeaf(EIG_MAT* serialize, int* T_ORG, double* double_buff){

    int step = 0;
    for (int i = 0; i < serialize->n_non_leaf; i++)
    {
        for (int m = 0; m < 5; m++)
        {

            memcpy(double_buff + step, serialize->Q0_nonleaf[i]->QC[m], sizeof(double) * serialize->Q0_nonleaf[i]->qcSizes[m].first * serialize->Q0_nonleaf[i]->qcSizes[m].second);
            step += serialize->Q0_nonleaf[i]->qcSizes[m].first * serialize->Q0_nonleaf[i]->qcSizes[m].second;
        }

        memcpy(double_buff + step, serialize->Q0_nonleaf[i]->J, sizeof(double) * serialize->Q0_nonleaf[i]->JSize.first * serialize->Q0_nonleaf[i]->JSize.second);
        step += serialize->Q0_nonleaf[i]->JSize.first * serialize->Q0_nonleaf[i]->JSize.second;

        memcpy(double_buff + step, serialize->Q0_nonleaf[i]->G, sizeof(double) * serialize->Q0_nonleaf[i]->GSize.first * serialize->Q0_nonleaf[i]->GSize.second);
        step += serialize->Q0_nonleaf[i]->GSize.first * serialize->Q0_nonleaf[i]->GSize.second;

        memcpy(double_buff + step, serialize->Q0_nonleaf[i]->I, sizeof(double) * serialize->Q0_nonleaf[i]->ISize.first * serialize->Q0_nonleaf[i]->ISize.second);
        step += serialize->Q0_nonleaf[i]->ISize.first * serialize->Q0_nonleaf[i]->ISize.second;

        memcpy(double_buff + step, serialize->Q0_nonleaf[i]->v2c, sizeof(double) * serialize->Q0_nonleaf[i]->v2cSize.first * serialize->Q0_nonleaf[i]->v2cSize.second);
        step += serialize->Q0_nonleaf[i]->v2cSize.first * serialize->Q0_nonleaf[i]->v2cSize.second;
    }

    step = 0;
    for (int i = 0; i < serialize->n_non_leaf; i++)
    {

        memcpy(T_ORG + step, serialize->Q0_nonleaf[i]->T, sizeof(int) * serialize->Q0_nonleaf[i]->TSize.first * serialize->Q0_nonleaf[i]->TSize.second);
        step += serialize->Q0_nonleaf[i]->TSize.first * serialize->Q0_nonleaf[i]->TSize.second;

        memcpy(T_ORG + step, serialize->Q0_nonleaf[i]->Org, sizeof(int) * serialize->Q0_nonleaf[i]->qcSizes[5].first * serialize->Q0_nonleaf[i]->qcSizes[5].second);
        step += serialize->Q0_nonleaf[i]->qcSizes[5].first * serialize->Q0_nonleaf[i]->qcSizes[5].second;
    }
    
}

void deserialzeQNonLeaf(EIG_MAT* deserialze, int* T_ORG, double* double_buff){

    int step=0;
    for (int i = 0; i < deserialze->n_non_leaf; i++)
    {
        for (int m = 0; m < 5; m++)
        {   
            deserialze->Q0_nonleaf[i]->QC[m] = (double_buff+step);
            step += deserialze->Q0_nonleaf[i]->qcSizes[m].first*deserialze->Q0_nonleaf[i]->qcSizes[m].second;
        }
        
        deserialze->Q0_nonleaf[i]->J = (double_buff+step);
        step+= deserialze->Q0_nonleaf[i]->JSize.first * deserialze->Q0_nonleaf[i]->JSize.second;

        deserialze->Q0_nonleaf[i]->G = (double_buff+step);
        step+= deserialze->Q0_nonleaf[i]->GSize.first * deserialze->Q0_nonleaf[i]->GSize.second;

        deserialze->Q0_nonleaf[i]->I = (double_buff+step);
        step+= deserialze->Q0_nonleaf[i]->ISize.first * deserialze->Q0_nonleaf[i]->ISize.second;

        deserialze->Q0_nonleaf[i]->v2c = (double_buff+step);
        step+= deserialze->Q0_nonleaf[i]->v2cSize.first * deserialze->Q0_nonleaf[i]->v2cSize.second;
    }

    step = 0;
    for (int i = 0; i < deserialze->n_non_leaf; i++)
    {
        deserialze->Q0_nonleaf[i]->T = (T_ORG + step);
        step+= deserialze->Q0_nonleaf[i]->TSize.first*deserialze->Q0_nonleaf[i]->TSize.second;

        deserialze->Q0_nonleaf[i]->Org = (T_ORG+step);
        step+= deserialze->Q0_nonleaf[i]->qcSizes[5].first*deserialze->Q0_nonleaf[i]->qcSizes[5].second;
    }
    
}


// double *send_buf = new double[q0_size_second * (12 + 10 + 4)];
// for (int j = 0; j < q0_size_second; j++)
// {
//     int buff_indx = 0;
//     buff_indx = j * (12 + 8 + 4);
//     int step = 0;
//     for (int k = 0; k < 6; k++)
//     {
//         send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->qcSizes[k].first;
//         step += 1;
//         send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->qcSizes[k].second;
//         step += 1;
//     }

//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->JSize.first;
//     step += 1;
//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->JSize.second;
//     step += 1;

//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->GSize.first;
//     step += 1;
//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->GSize.second;
//     step += 1;

//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->ISize.first;
//     step += 1;
//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->ISize.second;
//     step += 1;

//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->v2cSize.first;
//     step += 1;
//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->v2cSize.second;
//     step += 1;

//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->TSize.first;
//     step += 1;
//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->TSize.second;
//     step += 1;

//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->n;
//     step += 1;
//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->n1;
//     step += 1;
//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->n2;
//     step += 1;
//     send_buf[buff_indx + step] = Q0[indx]->Q0_nonleaf[j]->n3;
//     step += 1;
// }

// SEND_Error = MPI_Send(send_buf, q0_size_second * (12 + 10 + 4), MPI_DOUBLE, send_to, 0, process_grid);
// if (SEND_Error != MPI_SUCCESS)
// {
//     cout << "SendError at 670"
//          << "\n";
// }

// SEND_Error = MPI_Send(Lam[indx], LamSize, MPI_DOUBLE, send_to, 0, process_grid);
// if (SEND_Error != MPI_SUCCESS)
// {
//     cout << "SendError at 670"
//          << "\n";
// }
// delete[] send_buf;

// for (int l = 0; l < q0_size_second; l++)
// {
//     nonleaf *serialize = Q0[indx]->Q0_nonleaf[l];

//     int total_sz = 0;
//     for (int m = 0; m < 6; m++)
//     {
//         total_sz += serialize->qcSizes[m].first * serialize->qcSizes[m].second;
//     }

//     total_sz += ((serialize->JSize.first * serialize->JSize.second) + (serialize->GSize.first * serialize->GSize.second));
//     total_sz += ((serialize->ISize.first * serialize->ISize.second) + (serialize->v2cSize.first * serialize->v2cSize.second));
//     total_sz += ((serialize->TSize.first * serialize->TSize.second) + 4);

//     double *send_buf_doubles = new double[total_sz];
//     int step = 0;
//     for (int m = 0; m < 6; m++)
//     {
//         if (m != 5)
//         {
//             memcpy(send_buf_doubles + step, serialize->QC[m], sizeof(double) * serialize->qcSizes[m].first * serialize->qcSizes[m].second);
//         }
//         step += serialize->qcSizes[m].first * serialize->qcSizes[m].second;
//         assert(step < total_sz);
//     }

//     memcpy(send_buf_doubles + step, serialize->J, sizeof(double) * serialize->JSize.first * serialize->JSize.second);
//     step += serialize->JSize.first * serialize->JSize.second;
//     assert(step < total_sz);

//     memcpy(send_buf_doubles + step, serialize->G, sizeof(double) * serialize->GSize.first * serialize->GSize.second);
//     step += serialize->GSize.first * serialize->GSize.second;
//     assert(step < total_sz);

//     memcpy(send_buf_doubles + step, serialize->I, sizeof(double) * serialize->ISize.first * serialize->ISize.second);
//     step += serialize->ISize.first * serialize->ISize.second;
//     assert(step < total_sz);

//     memcpy(send_buf_doubles + step, serialize->v2c, sizeof(double) * serialize->v2cSize.first * serialize->v2cSize.second);
//     step += serialize->v2cSize.first * serialize->v2cSize.second;
//     assert(step < total_sz);

//     // // memcpy(send_buf_doubles+step, serialize->T, sizeof(double)*serialize->TSize.first*serialize->TSize.second);
//     step += serialize->TSize.first * serialize->TSize.second;
//     assert(step < total_sz);

//     SEND_Error = MPI_Send(send_buf_doubles, total_sz, MPI_DOUBLE, send_to, 0, process_grid);

//     if (SEND_Error != MPI_SUCCESS)
//     {
//         cout << "SendError at 670"
//              << "\n";
//     }

//     delete[] send_buf_doubles;

//     // send T and org seperately;
//     int *T_ORG = new int[serialize->qcSizes[5].first * serialize->qcSizes[5].second + serialize->TSize.first * serialize->TSize.second];
//     int offset = 0;

//     memcpy(T_ORG + offset, serialize->Org, sizeof(int) * serialize->qcSizes[5].first * serialize->qcSizes[5].second);
//     offset += serialize->qcSizes[5].first * serialize->qcSizes[5].second;

//     memcpy(T_ORG + offset, serialize->T, sizeof(int) * serialize->TSize.first * serialize->TSize.second);
//     offset += serialize->qcSizes[5].first * serialize->qcSizes[5].second;

//     SEND_Error = MPI_Send(T_ORG, serialize->qcSizes[5].first * serialize->qcSizes[5].second + serialize->TSize.first * serialize->TSize.second, MPI_INT, send_to, 0, process_grid);
//     if (SEND_Error != MPI_SUCCESS)
//     {
//         cout << "SendError at 670"
//              << "\n";
//     }
//     delete[] T_ORG;
// }