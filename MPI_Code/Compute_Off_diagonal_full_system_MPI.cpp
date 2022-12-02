//
// Created by phyzch on 6/30/20.
//
# include"../system.h"
#include "../util.h"
using namespace std;
double inter_detector_coupling_scaling=10;

void full_system::compute_offdiagonal_part_MPI(){
    int dmat_off_num, dmat_off_num_sum;
    int nonadiabatic_off_num, nonadiabatic_off_num_sum;

    vector<double> d_off_mat;
    vector<int> d_off_irow;
    vector<int> d_off_icol;
    compute_dmat_off_diagonal_matrix_in_full_matrix_MPI(d_off_mat,d_off_irow,d_off_icol);


    vector<double> nonadiabatic_off_mat;
    vector<int> nonadiabatic_off_irow;
    vector<int> nonadiabatic_off_icol;
    compute_nonadiabatic_offdiagonal_matrix_full_system(nonadiabatic_off_mat, nonadiabatic_off_irow, nonadiabatic_off_icol);

    // we have to rearrange off-diagonal_matrix in full_system to make sure irow is in  corresponding process.
    //Also we have to compute offnum, matnum
    combine_offdiagonal_term(d_off_mat,d_off_irow,d_off_icol,
                             nonadiabatic_off_mat, nonadiabatic_off_irow, nonadiabatic_off_icol);


    offnum= matnum-matsize;
    // compute total_matnum, total_offnum, matnum_each_process, offnum_each_process.
    mat_offnum_each_process= new int [num_proc];
    matnum_each_process= new int [num_proc];
    matnum_offset_each_process = new int [num_proc];
    MPI_Allgather(&matnum,1,MPI_INT,&matnum_each_process[0],1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&offnum,1,MPI_INT,&mat_offnum_each_process[0],1,MPI_INT,MPI_COMM_WORLD);
    int i;
    matnum_offset_each_process[0]=0;
    for(i=1;i<num_proc;i++){
        matnum_offset_each_process[i]= matnum_offset_each_process[i-1] + matnum_each_process[i-1];
    }

    MPI_Allreduce(&matnum,&total_matnum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&offnum,&total_offnum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    dmat_off_num = d_off_mat.size();
    MPI_Allreduce(&dmat_off_num, &dmat_off_num_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    nonadiabatic_off_num = nonadiabatic_off_mat.size();
    MPI_Allreduce(&nonadiabatic_off_num, &nonadiabatic_off_num_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(my_id == 0){
        output << "anharmonic coupling number in dimer : " << dmat_off_num_sum << endl;
        output << "nonadiabatic coupling number in dimer: " << nonadiabatic_off_num_sum << endl;
    }


}

void full_system:: combine_offdiagonal_term(
        vector<double> & d_off_mat, vector<int> & d_off_irow, vector<int> & d_off_icol,
        vector<double> & nonadiabatic_off_mat, vector<int> & nonadiabatic_off_irow, vector<int> & nonadiabatic_off_icol){
    // combine 3 part of off-diagonal term together and construct sdindex, sdmode, sdnum

    int m,i;
    matnum=matsize;
    int size;


    matnum = matnum + d_off_mat.size() + nonadiabatic_off_mat.size();  // off_mat is for off-diagonal coupling. d_d_mat is for coupling between detector.

    mat.reserve(matnum);
    irow.reserve(matnum);
    icol.reserve(matnum);
    int matindex = matsize;


    // push off-diagonal coupling terms into mat.
    size= d_off_mat.size();
    for(i=0;i<size;i++){
        mat.push_back(d_off_mat[i]);
        irow.push_back(d_off_irow[i]);
        icol.push_back(d_off_icol[i]);
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

void full_system:: rearrange_off_diagonal_term(vector < double > & sys_detector_mat,vector  <int> & sys_detector_irow,
                                               vector<int> & sys_detector_icol){
    /*
     rearrange term for off-diagonal part of mat, irow, icol. to make sure for each element, its irow is in
     corresponding process. Although called sys_detector_mat, it can be used for all off-diagonal part for full_matrix.
    */

    int i,j;
    int pc_index;
    int irow_index;
    int size;
    int index;
    int size_per_pc;
    size = sys_detector_irow.size();
    int * sys_detector_coup_to_send_pc_number = new int  [num_proc];  // number of element to send each pc
    int * sys_detector_coup_to_send_pc_displ= new int  [num_proc];  // displs for each pc.

    vector <int> * sys_detector_irow_to_send_each_pc = new vector <int> [num_proc];  // record irow to send for each pc
    vector <int> * sys_detector_icol_to_send_each_pc = new vector <int> [num_proc];  // record icol to send for each pc
    vector <double> * sys_detector_mat_to_send_each_pc = new vector <double> [num_proc]; // record mat to send for each pc

    int * sys_detector_irow_to_send = new int  [size];  //irow ordered by pc_index
    int * sys_detector_icol_to_send= new int  [size]; // icol ordered by pc_index
    double * sys_detector_mat_to_send= new double  [size];

    for(i=0;i<num_proc;i++){
        sys_detector_coup_to_send_pc_number[i] = 0;
    }
    for(i=0;i<size;i++){
        irow_index=sys_detector_irow[i];
        for(j=1;j<num_proc;j++){
            if(irow_index < matsize_offset_each_process[j]){
                break;
            }
        }
        pc_index= j-1;
        sys_detector_coup_to_send_pc_number[pc_index]++;
        sys_detector_irow_to_send_each_pc[pc_index].push_back(irow_index);
        sys_detector_icol_to_send_each_pc[pc_index].push_back(sys_detector_icol[i]);
        sys_detector_mat_to_send_each_pc[pc_index].push_back(sys_detector_mat[i]);
    }

    sys_detector_coup_to_send_pc_displ[0]=0;
    for(i=1;i<num_proc;i++){
        sys_detector_coup_to_send_pc_displ[i] = sys_detector_coup_to_send_pc_displ[i-1] +
                                                sys_detector_coup_to_send_pc_number[i-1];
    }

    index=0;
    for(i=0;i<num_proc;i++){
        size_per_pc = sys_detector_irow_to_send_each_pc[i].size();
        for(j=0;j<size_per_pc;j++){
            sys_detector_irow_to_send[index] = sys_detector_irow_to_send_each_pc[i][j] ;
            sys_detector_icol_to_send[index] = sys_detector_icol_to_send_each_pc[i][j];
            sys_detector_mat_to_send[index] = sys_detector_mat_to_send_each_pc[i][j];
            index++;
        }
    }

    int * sys_detector_coup_to_recv_num ;
    int * sys_detector_coup_to_recv_displ;
    int * sys_detector_coup_irow_to_recv;
    int * sys_detector_coup_icol_to_recv;
    double * sys_detector_coup_mat_to_recv;
    int  total_recv_num ;


    sys_detector_coup_to_recv_num = new int [num_proc];
    for(i=0;i<num_proc;i++){
        sys_detector_coup_to_recv_num[i]=0;
    }
    sys_detector_coup_to_recv_displ= new int [num_proc];
    MPI_Alltoall(sys_detector_coup_to_send_pc_number, 1,MPI_INT,
                 sys_detector_coup_to_recv_num, 1 , MPI_INT , MPI_COMM_WORLD);
    total_recv_num=0;
    for(i=0;i<num_proc;i++){
        total_recv_num = total_recv_num + sys_detector_coup_to_recv_num[i];
    }
    sys_detector_coup_to_recv_displ[0]=0;
    for(i=1;i<num_proc;i++){
        sys_detector_coup_to_recv_displ[i]= sys_detector_coup_to_recv_displ[i-1] + sys_detector_coup_to_recv_num[i-1];
    }
    sys_detector_coup_irow_to_recv = new int [total_recv_num];
    sys_detector_coup_icol_to_recv = new int [total_recv_num];
    sys_detector_coup_mat_to_recv = new double [total_recv_num];

    //  send irow, icol, mat
    MPI_Alltoallv(&sys_detector_irow_to_send[0],sys_detector_coup_to_send_pc_number,
                  sys_detector_coup_to_send_pc_displ,MPI_INT,&sys_detector_coup_irow_to_recv[0],
                  sys_detector_coup_to_recv_num, sys_detector_coup_to_recv_displ,MPI_INT,MPI_COMM_WORLD);

    MPI_Alltoallv(&sys_detector_icol_to_send[0],sys_detector_coup_to_send_pc_number,
                  sys_detector_coup_to_send_pc_displ,MPI_INT,&sys_detector_coup_icol_to_recv[0],
                  sys_detector_coup_to_recv_num,sys_detector_coup_to_recv_displ,MPI_INT,MPI_COMM_WORLD);

    MPI_Alltoallv(&sys_detector_mat_to_send[0],sys_detector_coup_to_send_pc_number,
                  sys_detector_coup_to_send_pc_displ, MPI_DOUBLE, &sys_detector_coup_mat_to_recv[0],
                  sys_detector_coup_to_recv_num,sys_detector_coup_to_recv_displ,MPI_DOUBLE,MPI_COMM_WORLD);

    // assign new value to mat, irow, icol
    sys_detector_mat.clear();
    sys_detector_mat.reserve(total_recv_num);
    sys_detector_irow.clear();
    sys_detector_irow.reserve(total_recv_num);
    sys_detector_icol.clear();
    sys_detector_icol.reserve(total_recv_num);

    for(i=0;i<total_recv_num;i++){
        sys_detector_mat.push_back(sys_detector_coup_mat_to_recv[i]);
        sys_detector_irow.push_back(sys_detector_coup_irow_to_recv[i]);
        sys_detector_icol.push_back(sys_detector_coup_icol_to_recv[i]);
    }


    // free the space
    delete [] sys_detector_coup_to_send_pc_number;
    delete [] sys_detector_coup_to_send_pc_displ;

    delete [] sys_detector_irow_to_send_each_pc;
    delete [] sys_detector_icol_to_send_each_pc;
    delete [] sys_detector_mat_to_send_each_pc;

    delete [] sys_detector_irow_to_send;
    delete [] sys_detector_icol_to_send;
    delete [] sys_detector_mat_to_send;


    delete [] sys_detector_coup_to_recv_num;
    delete [] sys_detector_coup_to_recv_displ;

    delete [] sys_detector_coup_irow_to_recv;
    delete [] sys_detector_coup_icol_to_recv;
    delete [] sys_detector_coup_mat_to_recv;

}

// for dmat off diagonal coupling.
void full_system::compute_dmat_off_diagonal_matrix_in_full_matrix_MPI(vector < double > & d_off_mat,vector  <int> & d_off_irow, vector<int> & d_off_icol){
    for(int i=0;i<s.electronic_state_num; i++){
        compute_dmat_off_diagonal_matrix_in_full_matrix_one_dmat_MPI(i,d_off_mat,d_off_irow,d_off_icol);
    }
    rearrange_off_diagonal_term(d_off_mat,d_off_irow,d_off_icol);
}

void full_system::compute_dmat_off_diagonal_matrix_in_full_matrix_one_dmat_MPI(int index,vector < double > & d_off_mat,vector  <int> & d_off_irow, vector<int> & d_off_icol){
    // we just use q_index easily compute the result.
    vector<quotient_state> * dlist_ptr;
    quotient_state * d_ptr;
    // index is index for monomer. We include anharmonic coupling in each monomer.
    if(index ==0){
        dlist_ptr = &(d1list);
    }
    else{
        dlist_ptr = &(d2list);
    }

    int i,j,k,l,p;
    int m,n;
    int size= (*dlist_ptr).size();
    int q_index_list_size;
    for(m=0;m<size;m++){
        d_ptr =&((*dlist_ptr)[m]);
        q_index_list_size = (*d_ptr).q_index_list.size();
        for(n=0;n<q_index_list_size;n++){
            i= (*d_ptr).q_index_list[n][0];   // index in dmat  (diag)
            j = (*d_ptr).q_index_list[n][1];  // index in dmat  (diag)
            k=  (*d_ptr).q_index_list[n][2];  // index in mat (diagonal)
            l= (*d_ptr).q_index_list[n][3];   // index in mat (diagonal)
            p= (*d_ptr).q_index_list[n][4];   // index in dmat (off_diag)  p in dmat link i,j. or dirow[p] = i, dicol[p]=j

            if(i!=j){
                // this is off-diagonal element
                d_off_mat.push_back((*d_ptr).dmat_value_list[n]);
                d_off_irow.push_back(k);
                d_off_icol.push_back(l);

                // we have to count a_{ji} to enable parallel version of SUR algorithm
                d_off_mat.push_back( (*d_ptr).dmat_value_list[n] );
                d_off_irow.push_back(l);
                d_off_icol.push_back(k);
            }
        }
    }
}


