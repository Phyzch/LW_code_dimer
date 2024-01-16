//
// Created by phyzch on 6/29/20.
//
#include "../system.h"
#include "../util.h"

using namespace std;
void full_system::construct_dimer_Hamiltonian_matrix_with_energy_window_MPI() {
    int i;

    // compute diagonal part of dimer Hamiltonian.
    compute_sstate_dstate_diagpart_dirow_dicol_MPI();

    // sort (exciton_state, vib_mode1, vib_mode2) into groups representing states with the same (exciton_state, vib_mode1) or (exciton_state, vib_mode2).
    construct_quotient_state_all_MPI();

    compute_full_Hamiltonian_offdiagonal_part_MPI();
}

void full_system:: compute_sstate_dstate_diagpart_dirow_dicol_MPI(){
   // compute diagonal part and exciton_state_index_list, vibrational_state_index_list for full_system.
    int i,j,k,l;
    int mat_index = 0;
    double energy;
    // each process responsible to couple monomer 1's state with monomer 2's state to form state in x,y.
    int vsize_d1 = d.total_monomer_mat_size[0] / num_proc;
    int vsize_d2 = d.total_monomer_mat_size[1] / num_proc;

    // we divide the whole matrix into different process according to the quantum number in monomer 1.
    int begin_index_d1 = my_id * vsize_d1;
    int end_index_d1;
    if(my_id != num_proc - 1){
        end_index_d1= (my_id +1) * vsize_d1;
    }
    else{
        end_index_d1 = d.total_monomer_mat_size[0];
    }

    int initial_state_d1_global_index = vsize_d1 * d.initial_state_pc_id[0] + d.initial_state_index[0];
    int initial_state_d2_global_index = vsize_d2 * d.initial_state_pc_id[1] + d.initial_state_index[1];

    initial_dimer_state_pc_id = d.initial_state_pc_id[0]; // states are sort according to d1 state.

    // exciton_state_index_list : dimer vibrational states' exciton state ( denote different electronic state or potential energy surface)
    // vibrational_state_index_list : dimer vibrational states' vibrational number in each monomer (vibrational_state_index_list[0], vibrational_state_index_list[1])
    vibrational_state_index_list = new vector <int> [s.exciton_state_num];

    for(j= begin_index_d1 ; j< end_index_d1;j++){ // index in monomer 1 (monomer 1)
        for (k=0;k<d.total_monomer_mat_size[1]; k++) {  // index in monomer 2 (monomer 2)

            for(i=0;i<s.tlmatsize;i++){ // index for different PES.

                if( j == initial_state_d1_global_index and k== initial_state_d2_global_index and i==0 ){
                    if (my_id != initial_dimer_state_pc_id){
                        printf("wrong. dimer pc id should be equal to detector1 state id.");
                        MPI_Abort(MPI_COMM_WORLD, -25);
                    }
                    initial_dimer_state_index = mat_index;
                }

                energy = s.electronic_state_energy[i] + d.dmat_diagonal_global0[j] + d.dmat_diagonal_global1[k]; // energy of electronic state (s.tlmat[i]) + energy in two monomer(d.monomer1_vib_state_energy_all_pc[j] + d.monomer2_vib_state_energy_all_pc[k]).
                exciton_state_index_list.push_back(i);
                vibrational_state_index_list[0].push_back(j); // vibrational_state_index_list record monomer index across process (global index , not index in each process)
                vibrational_state_index_list[1].push_back(k);

                mat.push_back(energy);  // diagonal part is energy of vibrational state
                irow.push_back(mat_index); //mat_index is local, we have to re-compute it after all process compute its own irow, icol. See code below
                icol.push_back(mat_index);
                mat_index = mat_index + 1;

            }
        }
    }

    // communicate about initial dimer state index.  dimer_state_pc_id is the same in all processes.
    MPI_Bcast(&initial_dimer_state_index, 1, MPI_INT, initial_dimer_state_pc_id, MPI_COMM_WORLD);

    // construct exciton_state_index_list_all, vibrational_state_index_list_all; for future convenience.
    vibrational_state_index_list_all = new vector<int> [s.exciton_state_num];
    for(j= 0 ; j< d.total_monomer_mat_size[0]; j++){ // index in monomer 1 (monomer 1)
        for (k=0;k<d.total_monomer_mat_size[1]; k++) {  // index in monomer 2 (monomer 2)
            for(i=0;i<s.tlmatsize;i++){

                exciton_state_index_list_all.push_back(i);
                vibrational_state_index_list_all[0].push_back(j); // vibrational_state_index_list record  global index of element in monomer
                vibrational_state_index_list_all[1].push_back(k);

            }
        }
    }


    // you have to re-assign the irow, icol in each process for global matrix.:
    matsize = mat.size();
    matsize_each_process = new int [num_proc];
    matsize_offset_each_process = new int [num_proc];
    MPI_Allgather(&matsize, 1, MPI_INT, &matsize_each_process[0], 1, MPI_INT,MPI_COMM_WORLD);

    matsize_offset_each_process[0] = 0;
    for (i = 1; i < num_proc; i++) {
        matsize_offset_each_process[i] = matsize_offset_each_process[i - 1] + matsize_each_process[i-1];
    }
    total_matsize=0;
    for(i=0;i<num_proc;i++){
        total_matsize= total_matsize + matsize_each_process[i];
    }

    int offset;
    offset = matsize_offset_each_process[my_id];
    // each process re-assign irow, icol according to offset. Now irow, icol is index in global.
    for(i=0;i<matsize;i++){
        irow[i]= irow[i] + offset;
        icol[i]= icol[i] + offset;
    }

}


// initialize the wavefunction of our system.
// This code uses function  etot(). MPI_version of etot() is not finished yet.
void full_system::Initial_state_MPI() {
    // Construct initial wavefunction for full matrix
    // we should use information stored in monomer_qn_list0, monomer_qn_list1.
    double x1, y1, x2, y2, x3, y3;
    int m, i;
    norm = 0;
    x.reserve(matsize);
    y.reserve(matsize);

    double value;
    for (i = 0; i < matsize; i++) {
        // (x1 + i y1) (x2 + i y2) (x3 + i y3)
        x1 = s.x_electronic[exciton_state_index_list[i]]; // here we only have two choice: state initially in electronic state 1 or in electronic state 2.
        y1 = s.y_electronic[exciton_state_index_list[i]];
        x2 = d.xd_all[0][vibrational_state_index_list[0][i]];
        y2 = d.yd_all[0][vibrational_state_index_list[0][i]];
        x3 = d.xd_all[1][vibrational_state_index_list[1][i]];
        y3 = d.yd_all[1][vibrational_state_index_list[1][i]];
        value = x1*(x2*x3 - y2*y3) - y1*(x2*y3 + x3*y2);
        x.push_back(value);
        value = y1*(x2*x3 - y2*y3) + x1*(x2*y3 + x3*y2);
        y.push_back(value);
        norm = norm + pow(x[i], 2) + pow(y[i], 2);
        if(isnan(norm)){
            cout<< "Norm is not a number" << endl;
        }
    }

    MPI_Allreduce(&norm, &total_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    total_norm = 1 / sqrt(total_norm);
    for (i = 0; i < matsize; i++) {
        x[i] = x[i] * total_norm;
        y[i] = y[i] * total_norm;
    }

}
void full_system::shift_mat(){
    int m;
    int i;
    // MPI_version of computing energy is not finished.
    double *hx, *hy;
    hx = new double[matsize];  // we have to rewrite this part. hx should be allocated space outside the function.
    hy = new double[matsize];
    etot_MPI(hx,hy);  // calculate total energy of system.  MPI_version is not finished yet.

    if(my_id==0) {
        output << "Total energy before SUR algorithm shifting  " << total_energy / cf << endl;
        output << "Intra-monomer coupling strength V set as: " << d.V_intra << endl;
        output << "Intra-monomer coupling scaling parameter a set as  " << d.a_intra << endl;
        output << "Initial Detector state is : " << endl;
        for (m = 0; m < s.exciton_state_num; m++) {
            for (i = 0; i < d.nmodes[m]; i++) {
                output << d.initial_vibrational_state[m][i] << " ";
            }
            output << endl;
        }
    }
    delete [] hx;
    delete [] hy;
}

