//
// Created by phyzch on 6/30/20.
//
# include"../system.h"
#include "../util.h"
using namespace std;
double inter_detector_coupling_scaling=10;

void full_system::compute_full_Hamiltonian_offdiagonal_part_MPI(){
    int i;
    int anharmonic_coupling_num, anharmonic_coupling_num_sum;
    int nonadiabatic_off_num, nonadiabatic_off_num_sum;

    // off-diagonal elements in Hamiltonian, due to anharmonic coupling between states in the same monomer
    vector<double> anharmonic_coupling_mat;
    vector<int> anharmonic_coupling_irow;
    vector<int> anharmonic_coupling_icol;
    compute_monomer_anharmonic_coupling_in_full_matrix_MPI(anharmonic_coupling_mat, anharmonic_coupling_irow,
                                                           anharmonic_coupling_icol);

    // off-diagonal elements in Hamiltonian, due to nonadiabatic coupling between states in different exciton states (potential energy surface)
    vector<double> nonadiabatic_off_mat;
    vector<int> nonadiabatic_off_irow;
    vector<int> nonadiabatic_off_icol;
    compute_nonadiabatic_offdiagonal_matrix_full_system(nonadiabatic_off_mat, nonadiabatic_off_irow, nonadiabatic_off_icol);

    //we have to rearrange off-diagonal_matrix in full_system to make sure irow is in  corresponding process.
    //Also we have to compute offnum, matnum
    combine_offdiagonal_term(anharmonic_coupling_mat, anharmonic_coupling_irow, anharmonic_coupling_icol,
                             nonadiabatic_off_mat, nonadiabatic_off_irow, nonadiabatic_off_icol);


    offnum = matnum-matsize;

    // compute total_matnum, total_offnum, matnum_each_process, offnum_each_process.
    mat_offnum_each_process = new int [num_proc];  // off-diagonal Hamiltonian elements number in each process
    matnum_each_process = new int [num_proc];  // Hamiltonian matrix number in each process
    matnum_offset_each_process = new int [num_proc]; // offset for matrix number in each process. Used for MPI communication.
    MPI_Allgather(&matnum, 1, MPI_INT, &matnum_each_process[0], 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&offnum, 1, MPI_INT, &mat_offnum_each_process[0], 1, MPI_INT, MPI_COMM_WORLD);

    matnum_offset_each_process[0] = 0;
    for(i=1;i<num_proc;i++){
        matnum_offset_each_process[i]= matnum_offset_each_process[i-1] + matnum_each_process[i-1];
    }

    MPI_Allreduce(&matnum,&total_matnum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&offnum,&total_offnum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    anharmonic_coupling_num = anharmonic_coupling_mat.size();
    MPI_Allreduce(&anharmonic_coupling_num, &anharmonic_coupling_num_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    nonadiabatic_off_num = nonadiabatic_off_mat.size();
    MPI_Allreduce(&nonadiabatic_off_num, &nonadiabatic_off_num_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(my_id == 0){
        output << "anharmonic coupling number in dimer : " << anharmonic_coupling_num_sum << endl;
        output << "nonadiabatic coupling number in dimer: " << nonadiabatic_off_num_sum << endl;
    }


}

void full_system:: combine_offdiagonal_term(
        vector<double> & anharmonic_coupling_mat, vector<int> & anharmonic_coupling_irow, vector<int> & anharmonic_coupling_icol,
        vector<double> & nonadiabatic_off_mat, vector<int> & nonadiabatic_off_irow, vector<int> & nonadiabatic_off_icol){
    // combine 3 part of off-diagonal term together and construct sdindex, sdmode, sdnum

    int m,i;
    int size;

    matnum = matsize;

    matnum = matnum + anharmonic_coupling_mat.size() + nonadiabatic_off_mat.size();  // off_mat is for off-diagonal coupling. d_d_mat is for coupling between monomer.

    mat.reserve(matnum);
    irow.reserve(matnum);
    icol.reserve(matnum);
    int matindex = matsize;


    // push off-diagonal coupling terms into mat.
    size = anharmonic_coupling_mat.size();
    for(i=0;i<size;i++){
        mat.push_back(anharmonic_coupling_mat[i]);
        irow.push_back(anharmonic_coupling_irow[i]);
        icol.push_back(anharmonic_coupling_icol[i]);
        matindex++;
    }

    size = nonadiabatic_off_mat.size();
    for(i=0;i<size;i++){
        mat.push_back(nonadiabatic_off_mat[i]);
        irow.push_back(nonadiabatic_off_irow[i]);
        icol.push_back(nonadiabatic_off_icol[i]);
        matindex ++;
    }

}

void full_system:: rearrange_matrix_element_in_different_pc(vector < double > & mat, vector  <int> & irow,
                                                            vector<int> & icol){
    /*
     rearrange term for off-diagonal part of mat, irow, icol. to make sure for each element, its irow is in
     corresponding process. Although called mat, it can be used for all off-diagonal part for full_matrix.
    */

    int i,j;
    int pc_index;
    int irow_index;
    int size;
    int index;
    int size_per_pc;

    size = irow.size();

    int * mat_element_number_to_send_pc = new int  [num_proc];  // number of element to send each pc
    int * mat_element_displacement_to_send_pc= new int  [num_proc];  // displacement for each pc.

    vector <int> * irow_to_send_diff_pc = new vector <int> [num_proc];  // record irow to send for each pc
    vector <int> * icol_to_send_diff_pc = new vector <int> [num_proc];  // record icol to send for each pc
    vector <double> * mat_value_to_send_diff_pc = new vector <double> [num_proc]; // record mat to send for each pc

    int * irow_to_send = new int  [size];  //irow ordered by pc_index
    int * icol_to_send = new int  [size]; // icol ordered by pc_index
    double * mat_value_to_send = new double  [size];

    for(i=0;i<num_proc;i++){
        mat_element_number_to_send_pc[i] = 0;
    }
    for(i=0;i<size;i++){
        irow_index = irow[i];

        // figure out process index.
        for(j = 1;j < num_proc; j++){
            if(irow_index < matsize_offset_each_process[j]){
                break;
            }
        }

        pc_index= j-1;

        mat_element_number_to_send_pc[pc_index]++;
        irow_to_send_diff_pc[pc_index].push_back(irow_index);
        icol_to_send_diff_pc[pc_index].push_back(icol[i]);
        mat_value_to_send_diff_pc[pc_index].push_back(mat[i]);
    }

    mat_element_displacement_to_send_pc[0] = 0;
    for(i=1;i<num_proc;i++){
        mat_element_displacement_to_send_pc[i] = mat_element_displacement_to_send_pc[i - 1] +
                                                 mat_element_number_to_send_pc[i - 1];
    }

    // sort irow, icol, mat_value according to the process index it is going to send to.
    index = 0;
    for(pc_index = 0; pc_index < num_proc; pc_index++){
        size_per_pc = irow_to_send_diff_pc[pc_index].size();
        for(j = 0; j < size_per_pc; j++){
            irow_to_send[index] = irow_to_send_diff_pc[pc_index][j] ;
            icol_to_send[index] = icol_to_send_diff_pc[pc_index][j];
            mat_value_to_send[index] = mat_value_to_send_diff_pc[pc_index][j];
            index++;
        }
    }

    int * matrix_element_number_to_recv ;
    int * matrix_element_to_recv_displacement;
    int * matrix_irow_to_recv;
    int * matrix_icol_to_recv;
    double * matrix_element_value_to_recv;
    int  total_recv_num ;


    matrix_element_number_to_recv = new int [num_proc];
    for(i = 0;i < num_proc; i++){
        matrix_element_number_to_recv[i]=0;
    }
    matrix_element_to_recv_displacement= new int [num_proc];

    MPI_Alltoall(mat_element_number_to_send_pc, 1, MPI_INT,
                 matrix_element_number_to_recv, 1 , MPI_INT , MPI_COMM_WORLD);

    total_recv_num=0;
    for(i = 0;i < num_proc; i++){
        total_recv_num = total_recv_num + matrix_element_number_to_recv[i];
    }

    matrix_element_to_recv_displacement[0] = 0;
    for(i = 1;i < num_proc; i++){
        matrix_element_to_recv_displacement[i] = matrix_element_to_recv_displacement[i - 1] + matrix_element_number_to_recv[i - 1];
    }

    matrix_irow_to_recv = new int [total_recv_num];
    matrix_icol_to_recv = new int [total_recv_num];
    matrix_element_value_to_recv = new double [total_recv_num];

    //  send irow, icol, mat
    MPI_Alltoallv(&irow_to_send[0], mat_element_number_to_send_pc,
                  mat_element_displacement_to_send_pc, MPI_INT, &matrix_irow_to_recv[0],
                  matrix_element_number_to_recv, matrix_element_to_recv_displacement, MPI_INT, MPI_COMM_WORLD);

    MPI_Alltoallv(&icol_to_send[0], mat_element_number_to_send_pc,
                  mat_element_displacement_to_send_pc, MPI_INT, &matrix_icol_to_recv[0],
                  matrix_element_number_to_recv, matrix_element_to_recv_displacement, MPI_INT, MPI_COMM_WORLD);

    MPI_Alltoallv(&mat_value_to_send[0], mat_element_number_to_send_pc,
                  mat_element_displacement_to_send_pc, MPI_DOUBLE, &matrix_element_value_to_recv[0],
                  matrix_element_number_to_recv, matrix_element_to_recv_displacement, MPI_DOUBLE, MPI_COMM_WORLD);

    // assign new value to mat, irow, icol
    mat.clear();
    mat.reserve(total_recv_num);
    irow.clear();
    irow.reserve(total_recv_num);
    icol.clear();
    icol.reserve(total_recv_num);

    for(i=0;i<total_recv_num;i++){
        mat.push_back(matrix_element_value_to_recv[i]);
        irow.push_back(matrix_irow_to_recv[i]);
        icol.push_back(matrix_icol_to_recv[i]);
    }


    // free the space
    delete [] mat_element_number_to_send_pc;
    delete [] mat_element_displacement_to_send_pc;

    delete [] irow_to_send_diff_pc;
    delete [] icol_to_send_diff_pc;
    delete [] mat_value_to_send_diff_pc;

    delete [] irow_to_send;
    delete [] icol_to_send;
    delete [] mat_value_to_send;


    delete [] matrix_element_number_to_recv;
    delete [] matrix_element_to_recv_displacement;

    delete [] matrix_irow_to_recv;
    delete [] matrix_icol_to_recv;
    delete [] matrix_element_value_to_recv;

}

// for anharmonic coupling in the same monomer coupling.
void full_system::compute_monomer_anharmonic_coupling_in_full_matrix_MPI(vector < double > & anharmonic_coupling_mat, vector  <int> & anharmonic_coupling_irow, vector<int> & anharmonic_coupling_icol){
    // compute anharmonic coupling within monomer in each exciton state.
    int i;
    for(i = 0; i < monomer_number; i++){
        compute_anharmonic_coupling_in_full_matrix_in_one_monomer_MPI(i, anharmonic_coupling_mat,
                                                                      anharmonic_coupling_irow,
                                                                      anharmonic_coupling_icol);
    }

    rearrange_matrix_element_in_different_pc(anharmonic_coupling_mat, anharmonic_coupling_irow,
                                             anharmonic_coupling_icol);
}

void full_system::compute_anharmonic_coupling_in_full_matrix_in_one_monomer_MPI(int monomer_index, vector < double > & anharmonic_coupling_mat, vector  <int> & anharmonic_coupling_irow, vector<int> & anharmonic_coupling_icol ){
    // we just use q_index to easily compute the result.
    vector<quotient_state> * monomer_state_list_ptr;
    quotient_state * monomer_state_ptr;
    //  We include anharmonic coupling in each monomer.
    if(monomer_index == 0){
        monomer_state_list_ptr = &(monomer1_quotient_state_list);
    }
    else{
        monomer_state_list_ptr = &(monomer2_quotient_state_list);
    }

    int i,j,k,l,p;
    int m,n;
    int monomer_state_list_size = (*monomer_state_list_ptr).size();
    int anharmonic_coupling_info_index_list_size;

    for(m = 0; m < monomer_state_list_size; m++){
        monomer_state_ptr = &((*monomer_state_list_ptr)[m]);
        anharmonic_coupling_info_index_list_size = (*monomer_state_ptr).anharmonic_coupling_info_index_list.size();

        for(n = 0; n < anharmonic_coupling_info_index_list_size; n++){
            i = (*monomer_state_ptr).anharmonic_coupling_info_index_list[n][0];   // monomer_index in vibrational_state_index_list (monomer1)
            j = (*monomer_state_ptr).anharmonic_coupling_info_index_list[n][1];  // monomer_index in vibrational_state_index_list (monomer2)
            k =  (*monomer_state_ptr).anharmonic_coupling_info_index_list[n][2];  // states in full system's basis set which corresponds to monomer state i.
            l = (*monomer_state_ptr).anharmonic_coupling_info_index_list[n][3];   // states in full system's basis set which corresponds to monomer state j.
            p = (*monomer_state_ptr).anharmonic_coupling_info_index_list[n][4];   // monomer_index in monomer_mat (off_diag)  p in monomer Hamiltonian link i,j. (monomer_irow[p] = i, monomer_icol[p] = j).

            if(i!=j){
                // this is off-diagonal element
                anharmonic_coupling_mat.push_back((*monomer_state_ptr).anharmonic_coupling_value_list[n]);
                anharmonic_coupling_irow.push_back(k);
                anharmonic_coupling_icol.push_back(l);

                // we have to count a_{ji} to enable parallel version of SUR algorithm
                anharmonic_coupling_mat.push_back((*monomer_state_ptr).anharmonic_coupling_value_list[n] );
                anharmonic_coupling_irow.push_back(l);
                anharmonic_coupling_icol.push_back(k);
            }
        }

    }
}


