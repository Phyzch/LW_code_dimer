//
// Created by phyzch on 10/28/22.
//
#include"system.h"
#include"util.h"
using namespace std;


void full_system:: compute_local_density_of_state(vector<vector<vector<int>>> & coupling_state_index_list,
                                                  vector<vector<vector<vector<int>>>> & coupling_state_qn_list,
                                                  vector<vector<vector<double>>> & coupling_state_strength_list,
                                                  vector<vector<vector<double>>> & coupling_state_energy_diff_list,
                                                  vector<vector<double>> & effective_coupling_number_list){

    d.compute_local_density_of_state(coupling_state_index_list, coupling_state_qn_list, coupling_state_strength_list,
                                   coupling_state_energy_diff_list, effective_coupling_number_list);
}

void detector:: compute_local_density_of_state(vector<vector<vector<int>>> & coupling_state_index_list,
                                               vector<vector<vector<vector<int>>>> & coupling_state_qn_list,
                                               vector<vector<vector<double>>> & coupling_state_strength_list,
                                               vector<vector<vector<double>>> & coupling_state_energy_diff_list,
                                               vector<vector<double>> & effective_coupling_number_list){
    int i,j,k;
    int monomer_index;

    vector<vector<int>> coupling_state_index;
    vector<vector<vector<int>>> coupling_state_qn;
    vector<vector<double>> coupling_state_strength;
    vector<vector<double>> coupling_state_energy_diff;
    vector<double> effective_coupling_number;
    for(monomer_index = 0; monomer_index < electronic_state_num; monomer_index ++ ){
        vector<int> coupling_state_index_each_monomer;
        vector<vector<int>> coupling_state_qn_each_monomer;
        vector<double> coupling_state_strength_each_monomer;
        vector<double> coupling_state_energy_diff_each_monomer;
        double effective_coupling_number_each_monomer;

        if (monomer_index == 0){
            compute_local_density_of_state_subroutine(monomer_index, dmat_diagonal_global0,
                                                      coupling_state_index_each_monomer, coupling_state_qn_each_monomer,
                                                      coupling_state_strength_each_monomer, coupling_state_energy_diff_each_monomer,
                                                      effective_coupling_number_each_monomer);
        }
        else{
            compute_local_density_of_state_subroutine(monomer_index, dmat_diagonal_global1,
                                                      coupling_state_index_each_monomer, coupling_state_qn_each_monomer,
                                                      coupling_state_strength_each_monomer, coupling_state_energy_diff_each_monomer,
                                                      effective_coupling_number_each_monomer);
        }

        coupling_state_index.push_back(coupling_state_index_each_monomer);
        coupling_state_qn.push_back(coupling_state_qn_each_monomer);
        coupling_state_strength.push_back(coupling_state_strength_each_monomer);
        coupling_state_energy_diff.push_back(coupling_state_energy_diff_each_monomer);
        effective_coupling_number.push_back(effective_coupling_number_each_monomer);
    }

    coupling_state_index_list.push_back(coupling_state_index);
    coupling_state_qn_list.push_back(coupling_state_qn);
    coupling_state_strength_list.push_back(coupling_state_strength);
    coupling_state_energy_diff_list.push_back(coupling_state_energy_diff);
    effective_coupling_number_list.push_back(effective_coupling_number);
}

void detector:: compute_local_density_of_state_subroutine(int monomer_index, vector<double> & state_energy_global_matrix,
                                                          vector<int> & coupling_state_index_list, vector<vector<int>> & coupling_state_qn_list,
                                                          vector<double> & coupling_state_strength, vector<double> & coupling_state_energy_diff,
                                                          double & effective_coupling_number)
                                                          {
    // state_index_in_global_matrix is index across all process.
    // We want to compute : (1) strength of coupling Vj  (2) vibrational state quantum number (v1, v2, v3) couple to given state
    // (3) energy difference between states coupled state j and state i  (4) effective coupling number L == 1/( 1 + (Eij / Vij)^2 )
    // same PES bool will decide which kind of coupling we want to compute : (1) coupling between diff PES if False (2) coupling within same PES if True
    // effective_coupling_number: Li = 1/ (1 + (Delta E/ V)^2)
    int i , j, k;

    int state_index_in_global_matrix;
    int state_proc_id;
    int state_index_in_proc;



    state_proc_id = initial_state_pc_id[monomer_index];
    state_index_in_proc = initial_state_index[monomer_index];
    state_index_in_global_matrix = dmatsize_offset_each_process[monomer_index][state_proc_id] + state_index_in_proc;


    double energy_diff;
    int coupling_state_index;
    vector<vector<int>> & vibration_qn = dv_all[monomer_index];

    // for record result
    int coupling_state_num;
    effective_coupling_number = 0;
    double effective_coupling_number_between_states;

    if (my_id == state_proc_id){
        // go through off diagonal coupling
        for( i = dmatsize[monomer_index] ; i<dmatnum[monomer_index]; i++){
            if(dirow[monomer_index][i] == state_index_in_global_matrix){
                energy_diff = abs(state_energy_global_matrix[state_index_in_global_matrix] - state_energy_global_matrix[dicol[monomer_index][i]] );

                // do not include state whose energy difference is small (polyad).
//                if( abs(energy_diff) <= 10 ){
//                    continue;
//                }
                effective_coupling_number_between_states = 1 / (1 + pow(energy_diff / dmat[monomer_index][i] , 2));
                effective_coupling_number = effective_coupling_number + effective_coupling_number_between_states ;
                coupling_state_index = dicol[monomer_index][i] ;

                if (effective_coupling_number_between_states > 0.01 ){
                    coupling_state_index_list.push_back(coupling_state_index );
                    coupling_state_energy_diff.push_back(energy_diff);
                    coupling_state_strength.push_back(dmat[monomer_index][i]);
                }

            }
        }

        coupling_state_num = coupling_state_index_list.size();
    }

    MPI_Bcast(&effective_coupling_number, 1, MPI_DOUBLE, state_proc_id, MPI_COMM_WORLD);
    broadcast_1d_vector<int>(coupling_state_index_list, coupling_state_num, MPI_INT, state_proc_id);
    broadcast_1d_vector<double>(coupling_state_strength, coupling_state_num, MPI_DOUBLE, state_proc_id);
    broadcast_1d_vector<double>(coupling_state_energy_diff, coupling_state_num, MPI_DOUBLE, state_proc_id);

    for(i=0; i < coupling_state_num; i++){
        coupling_state_index = coupling_state_index_list[i];
        coupling_state_qn_list.push_back(vibration_qn[coupling_state_index] );
    }

}
