//
// Created by phyzch on 7/19/20.
//
#include"../util.h"
#include"../system.h"
using namespace std;

// function is almost the same for the one Detector_Evolve_MPI
int full_system::construct_recvbuffer_index(){
    int i,j,k;
    int total_recv_count=0;
    int local_begin= matsize_offset_each_process[my_id];
    int local_end;
    if(my_id!=num_proc-1){
        local_end = matsize_offset_each_process[my_id+1];
    }
    else{
        local_end= total_matsize;
    }

    // --------- Allocate space for vector to receive -----------------------
    remoteVecCount = new int [num_proc];
    remoteVecPtr= new int [num_proc];
    remoteVecIndex = new int [offnum];
    for(i=0;i<num_proc;i++){
        remoteVecCount[i] = 0;
    }
    // ----------------------------------------------------------------


    vector<int> icol_copy = icol;
    sort(icol_copy.begin(),icol_copy.end());
    j=0;
    int prev_col=-1;
    int col_index;
    int col_pc_id;

    for(i = 0;i < matnum; i++){
        col_index = icol_copy[i];
        if(col_index > prev_col and (col_index < local_begin or col_index >= local_end) ){
            // find process id this icol belong to
            for(k=1;k<num_proc;k++){
                if(col_index < matsize_offset_each_process[k]){
                    break;
                }
            }
            col_pc_id= k-1;

            remoteVecCount[col_pc_id]++;
            remoteVecIndex[j]= col_index;

            j++;
        }
        prev_col= col_index;
    }

    remoteVecPtr[0]=0;
    for(i=1;i<num_proc;i++){
        remoteVecPtr[i]= remoteVecPtr[i-1] + remoteVecCount[i-1];
    }
    total_recv_count=0;
    for(i=0;i<num_proc;i++){
        total_recv_count = total_recv_count + remoteVecCount[i];
    }
    return total_recv_count;
}


void full_system:: prepare_evolution(){

    int m,i;
    int col_index_for_search;
    int * search_Ind;
    // You have to allocate space for tosendVecCount, tosendVecPtr.
    tosendVecCount = new int [num_proc];
    tosendVecPtr = new int [num_proc];

    to_recv_buffer_len = construct_recvbuffer_index();

    to_send_buffer_len = construct_send_buffer_index(remoteVecCount,remoteVecPtr,remoteVecIndex,
                                                    tosendVecCount,tosendVecPtr,tosendVecIndex);

    // we add extra space at the end of wave function array real_part_wave_func,imag_part_wave_func.
    // Which will store the wave function elements received from other process through MPI.
    real_part_wave_func.resize(matsize + to_recv_buffer_len);
    imag_part_wave_func.resize(matsize + to_recv_buffer_len);

    recv_real_wave_func = new double [to_recv_buffer_len];
    recv_imag_wave_func = new double [to_recv_buffer_len];
    send_real_wave_func = new double[to_send_buffer_len];
    send_imag_wave_func = new double [to_send_buffer_len];

    local_irow.reserve(matnum);  // row index for matrix multiplication in the single process.
    local_icol.reserve(matnum);  // col index for matrix multiplication in the single process.
    for(i = 0;i < matnum; i++){
        local_irow.push_back(irow[i] - matsize_offset_each_process[my_id]);
        col_index_for_search = icol[i];
        search_Ind = (int *)bsearch(&col_index_for_search,remoteVecIndex,to_recv_buffer_len,sizeof(int),compar);
        if (search_Ind != NULL){
            // this column element is in remote vector:
            local_icol.push_back(matsize + (search_Ind - remoteVecIndex));
        }
        else{
            // this column element is in local vector:
            local_icol.push_back(icol[i] - matsize_offset_each_process[my_id]);
        }
    }

}

void full_system:: full_system_SUR_one_step()
{
    int i;
    int irow_index, icol_index;
    update_x_y();
    // SUR_algorithm
    for(i=0;i<matnum;i++){
        irow_index = local_irow[i];
        icol_index= local_icol[i];
        real_part_wave_func[irow_index] = real_part_wave_func[irow_index] + mat[i] * imag_part_wave_func[icol_index];
    }
    for(i=0;i<matnum;i++){
        irow_index=local_irow[i];
        icol_index= local_icol[i];
        imag_part_wave_func[irow_index] = imag_part_wave_func[irow_index] - mat[i] * real_part_wave_func[icol_index];
    }
}

void full_system:: update_x_y(){
    int i;
    int local_index;
    // update tosend_xd for sending  to other process
    for(i=0;i<to_send_buffer_len;i++){
        local_index = tosendVecIndex[i] - matsize_offset_each_process[my_id];
        send_real_wave_func[i] = real_part_wave_func[local_index];
        send_imag_wave_func[i] = imag_part_wave_func[local_index];
    }

    MPI_Alltoallv(&send_real_wave_func[0], tosendVecCount, tosendVecPtr, MPI_DOUBLE,
                  &recv_real_wave_func[0], remoteVecCount, remoteVecPtr, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Alltoallv(&send_imag_wave_func[0], tosendVecCount, tosendVecPtr, MPI_DOUBLE,
                  &recv_imag_wave_func[0], remoteVecCount, remoteVecPtr, MPI_DOUBLE, MPI_COMM_WORLD);
    for(i=0;i<to_recv_buffer_len;i++){
        real_part_wave_func[matsize + i] = recv_real_wave_func[i];
        imag_part_wave_func[matsize + i] = recv_imag_wave_func[i];
    }
}

void full_system:: update_real_part(){
    int i;
    int local_index;
    for(i=0;i<to_send_buffer_len;i++){
        local_index = tosendVecIndex[i] - matsize_offset_each_process[my_id];
        send_real_wave_func [i] = real_part_wave_func [local_index];
    }

    MPI_Alltoallv(&send_real_wave_func[0], tosendVecCount, tosendVecPtr, MPI_DOUBLE,
                  &recv_real_wave_func[0], remoteVecCount, remoteVecPtr, MPI_DOUBLE, MPI_COMM_WORLD);
    for(i=0;i<to_recv_buffer_len;i++){
        real_part_wave_func[matsize + i] = recv_real_wave_func[i];
    }
}

void full_system::update_imag_part(){
    int i;
    int local_index;
    // update tosend_xd for sending  to other process
    for(i=0;i<to_send_buffer_len;i++){
        local_index = tosendVecIndex[i] - matsize_offset_each_process[my_id];
        send_imag_wave_func [i] = imag_part_wave_func [local_index];
    }

    MPI_Alltoallv(&send_imag_wave_func[0], tosendVecCount, tosendVecPtr, MPI_DOUBLE,
                  &recv_imag_wave_func[0], remoteVecCount, remoteVecPtr, MPI_DOUBLE, MPI_COMM_WORLD);

    for(i=0;i<to_recv_buffer_len;i++){
        imag_part_wave_func [matsize + i] = recv_imag_wave_func [i];
    }
}

void full_system:: evolve_wave_func_one_step(){
    // We can also use Chebyshev method here. But here we use SUR algorithm which gets the job done.
    int i, j, irow_index, icol_index;
    // update imaginary part of wave function received from other process.
    update_imag_part();

    // SUR algorithm. (https://doi.org/10.1016/0009-2614(94)01474-A)
    for(i=0;i<matnum;i++){
        irow_index = local_irow[i];
        icol_index= local_icol[i];
        real_part_wave_func[irow_index] = real_part_wave_func[irow_index] + mat[i] * imag_part_wave_func[icol_index];
    }

    // update real part of wave function received from other process.
    update_real_part();
    for(i=0;i<matnum;i++){
        irow_index=local_irow[i];
        icol_index= local_icol[i];
        imag_part_wave_func[irow_index] = imag_part_wave_func[irow_index] - mat[i] * real_part_wave_func[icol_index];
    }
}

void full_system::Normalize_wave_function(){
    int i;
    norm = 0;
    total_norm = 0;
    for(i=0;i<matsize;i++) {
        norm = norm + real_part_wave_func[i] * real_part_wave_func[i] + imag_part_wave_func[i] * imag_part_wave_func[i];
    }
    MPI_Allreduce(&norm,&total_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<matsize;i++){
        real_part_wave_func[i] = real_part_wave_func[i] / sqrt(total_norm);
        imag_part_wave_func[i] = imag_part_wave_func[i] / sqrt(total_norm);
    }
}