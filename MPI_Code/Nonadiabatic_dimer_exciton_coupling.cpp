//
// Created by phyzch on 11/30/22.
//

#include "../util.h"
#include "../system.h"
#include "../quotient_state.h"

int full_system::search_full_hamiltonian_state_index(int exciton_state, int monomer1_state, int monomer2_state){
    // exciton_state : index for electronic state
    // d_state: index for detector1
    // d_state: index for detector2

    // full_hamiltonian_state_index is index in full_matrix, which is the result we return. if full_hamiltonian_state_index = -1, this means the state is not found.

    int full_hamiltonian_state_index; // index in full matrix.

    vector <int> d1_mode;
    vector <int> d2_mode;

    d1_mode = d.monomer_vibrational_states_all[0][monomer1_state];
    d2_mode = d.monomer_vibrational_states_all[1][monomer2_state];

    int vsize_d1 = d.total_monomer_mat_size[0] / num_proc;

    // because full_system matrix are constructed with guarantee monomer1_state is sorted. Therefore,
    // we can use this attribute to find out in which proc the state locates.
    bool exist;
    int monomer_quotient_list_position;
    int full_hamiltonian_state_index_position;

    // search from d2_list_all.
    monomer_quotient_list_position = find_location_binarysearch_quotient_state(monomer2_quotient_state_list, exciton_state, d1_mode, exist);
    if(!exist){
        full_hamiltonian_state_index = -1; // no such state.
    }
    else{
        full_hamiltonian_state_index_position = binary_search_monomer_state_index(monomer2_quotient_state_list[monomer_quotient_list_position].monomer_state_index_list,
                                                                                  monomer2_state, exist);
        if(!exist){
            full_hamiltonian_state_index = -1;
        }
        else{
            full_hamiltonian_state_index = monomer2_quotient_state_list[monomer_quotient_list_position].full_hamiltonian_state_index_list[full_hamiltonian_state_index_position];

            // check if result is correct.
            if ((exciton_state_index_list_all[full_hamiltonian_state_index] != exciton_state) or (vibrational_state_index_list_all[0][full_hamiltonian_state_index] != monomer1_state) or (vibrational_state_index_list_all[1][full_hamiltonian_state_index] != monomer2_state ) ){
                printf("search state index in full system is wrong.");
                MPI_Abort(MPI_COMM_WORLD, -30);
            }

        }
    }

    return full_hamiltonian_state_index;
}

double compute_franck_condon_factor(double alpha, int m, int n){
    // compute franck_condon factor <m|alpha;n>. Here alpha is displacement operator.
    int nm_min = min(n,m);
    int l;
    double franck_condon_factor = 0;
    // prefacto = e^{-alpha^2/2} sqrt(n! m!) * alpha^{n+m}
    // std::tgamma(n+1) = n!
    double prefactor = exp(- pow(alpha,2)/2) * pow(alpha, n + m) * sqrt(std::tgamma(n+1)  * std::tgamma(m+1) );

    if (alpha!=0){
        for (l=0;l<=nm_min;l++){
            // sum: 1/(l! * (n-l)! * (m-l)!) * (-1)^{n-l} * alpha^{-2l}
            franck_condon_factor = franck_condon_factor + 1/( std::tgamma(l+1) * std::tgamma(n-l+1) * std::tgamma(m-l+1) ) * pow(-1, n-l) * pow(alpha, -2*l);
        }

        franck_condon_factor = franck_condon_factor * prefactor;
    }

    if (alpha ==0){
        if (m==n){
            franck_condon_factor = 1;
        }
        else{
            franck_condon_factor = 0;
        }
    }

    return franck_condon_factor;
}


void monomer::compute_franck_condon_factor_table() {
    // compute list of franck condon factor for nonadiabatic coupling.
    // different dof have different displacement (alpha).
    // results store in franck_condon_factor_table

    // max quantum number
    int max_qn = 0;
    int m;
    int i, j, k;
    double franck_condon_factor;
    double alpha;
    for(i=0;i<nmodes[0];i++){
        if (max_qn < nmax[0][i]) {
            max_qn = nmax[0][i];
        }
    }
    int monomer_number = 2;

    // include states from 0 to max_qn;
    max_qn = max_qn + 1;

    // compute franck condo factor for different monomers.
    franck_condon_factor_table = new double *** [monomer_number];
    for (m = 0; m < monomer_number; m++){
        franck_condon_factor_table[m] = new double **[nmodes[m]];
        for(j = 0;j < nmodes[m];j++){
            franck_condon_factor_table[m][j] = new double * [max_qn];
            for(k=0; k<max_qn; k++){
                franck_condon_factor_table[m][j][k] = new double [max_qn];
            }
        }
    }

    for (m=0; m < monomer_number ; m++ ){
        for(i=0;i<nmodes[m];i++){
            alpha = electron_phonon_coupling[m][i];
            for(j=0;j<max_qn;j++){
                for(k=0;k<max_qn;k++){
                    franck_condon_factor = compute_franck_condon_factor(alpha, j , k);
                    franck_condon_factor_table[m][i][j][k] = franck_condon_factor;
                }
            }
        }
    }


    // output franck condon factor table <i| alpha;j>
    if (my_id == 0){
        ofstream franck_condon_output;
        franck_condon_output.open(path + "franck_condon.txt");
        for(m=0;m<monomer_number;m++){
            for(i=0;i<nmodes[m];i++){
                franck_condon_output << "alpha: " << electron_phonon_coupling[m][i] << endl;
                for(j=0;j<max_qn;j++){
                    for(k=0;k<max_qn;k++){
                        franck_condon_output << franck_condon_factor_table[m][i][j][k] << " ";
                    }
                    franck_condon_output << endl;
                }
                franck_condon_output << endl;
            }

            franck_condon_output << endl;
        }

    }

}

void monomer::find_franck_condon_factor_for_monomer_states(){
    int i, j, k, l, m;
    double franck_condon_factor_prod;
    double franck_condon_factor;
    double energy_state_i, energy_state_j;
    double energy_difference;

    int state1_mode_qn, state2_mode_qn;

    int monmer_number = 2;

    // brute force method to find all monomer states
    for(m=0; m < monmer_number ; m++){
        vector<vector<int>> nonadiabatic_coupled_state_each_monomer;
        vector<vector<double>> nonadiabatic_coupled_state_each_monomer_fc;

        for(i = 0; i < total_monomer_mat_size[m] ; i++ ){
            vector<int> nonadiabatic_coupled_state_single_state;
            vector<double> nonadiabatic_coupled_state_single_state_fc;

            energy_state_i = 0;
            for(k = 0;k < nmodes[m]; k++){
                energy_state_i = energy_state_i + monomer_vibrational_states_all[m][i][k] * mfreq[m][k];
            }

            for(j = 0; j < total_monomer_mat_size[m]; j++ ){
                franck_condon_factor_prod = 1;

                for(k = 0; k < nmodes[m]; k++){
                    state1_mode_qn = monomer_vibrational_states_all[m][i][k];
                    state2_mode_qn = monomer_vibrational_states_all[m][j][k];

                    franck_condon_factor = franck_condon_factor_table[m][k][state1_mode_qn][state2_mode_qn];
                    franck_condon_factor_prod = franck_condon_factor_prod * franck_condon_factor;
                }

                energy_state_j = 0;
                for(k=0;k<nmodes[m];k++){
                    energy_state_j = energy_state_j + monomer_vibrational_states_all[m][j][k] * mfreq[m][k];
                }

                energy_difference = abs(energy_state_j - energy_state_i);

                if( abs(franck_condon_factor_prod) > Franck_condon_factor_cutoff or abs(franck_condon_factor_prod) > energy_difference * cutoff / nonadiabatic_coupling ){
                    nonadiabatic_coupled_state_single_state.push_back(j);
                    nonadiabatic_coupled_state_single_state_fc.push_back(franck_condon_factor_prod);
                }

            }
            nonadiabatic_coupled_state_each_monomer.push_back(nonadiabatic_coupled_state_single_state);
            nonadiabatic_coupled_state_each_monomer_fc.push_back(nonadiabatic_coupled_state_single_state_fc);
        }
        nonadiabatic_coupled_monomer_state.push_back(nonadiabatic_coupled_state_each_monomer);
        nonadiabatic_coupled_monomer_state_franck_condon.push_back(nonadiabatic_coupled_state_each_monomer_fc);

    }

}


void full_system:: compute_nonadiabatic_offdiagonal_matrix_full_system(vector < double > & nonadiabatic_off_mat, vector  <int> & nonadiabatic_off_irow, vector<int> & nonadiabatic_off_icol){
    int i,j,k;
    // exciton & monomer state index for full_sys state we are considering.
    int state_exciton_state_index;
    int state_monomer1_state_index, state_monomer2_state_index;

    // electronic & monomer state index for nonadiabatically coupled states
    int coupled_state_exciton_state_index;
    int coupled_state_monomer1_state_index, coupled_state_monomer2_state_index;

    int state_index_in_dimer;
    int coupled_state_index_in_dimer;

    // number of coupled vib states for given monomer
    int coupled_monomer1_state_num, coupled_monomer2_state_num;

    // value for off-diagonal coupling matrix element.
    double off_diagonal_matrix_ele;
    double franck_condon_factor_monomer1, franck_condon_factor_monomer2;

    double state_energy, coupled_state_energy;
    double state_energy_difference;

    bool include_nonadiabatic_coupling_bool; // bool to decide if including off diagonal coupling

    // compute franck condon factor <i| alpha;j>. result stored in d.franck_condon_factor_table
    d.compute_franck_condon_factor_table();

    // for each vib states in monomer, we find vib states coupled to them whose franck condon factor is larger than Franck_condon_factor_cutoff
    d.find_franck_condon_factor_for_monomer_states();

    for(i = 0; i < matsize; i++){
        state_index_in_dimer = irow[i];

        state_exciton_state_index = exciton_state_index_list[i];
        state_monomer1_state_index = vibrational_state_index_list[0][i]; // global index for monomer (monomer)
        state_monomer2_state_index = vibrational_state_index_list[1][i];

        state_energy = s.electronic_state_energy[state_exciton_state_index] + monomer1_vib_state_energy_all_pc[state_monomer1_state_index] + monomer2_vib_state_energy_all_pc[state_monomer2_state_index];

        if (state_exciton_state_index == 0){
            coupled_state_exciton_state_index = 1;
        }
        else{
            coupled_state_exciton_state_index = 0;
        }

        coupled_monomer1_state_num = d.nonadiabatic_coupled_monomer_state[0][state_monomer1_state_index].size();
        coupled_monomer2_state_num = d.nonadiabatic_coupled_monomer_state[1][state_monomer2_state_index].size();

        for(j = 0;j < coupled_monomer1_state_num;j++){
            coupled_state_monomer1_state_index = d.nonadiabatic_coupled_monomer_state[0][state_monomer1_state_index][j];
            for(k = 0; k < coupled_monomer2_state_num; k++){
                coupled_state_monomer2_state_index = d.nonadiabatic_coupled_monomer_state[1][state_monomer2_state_index][k];

                include_nonadiabatic_coupling_bool = false;

                // coupled state index in full matrix (dimer)
                // search full_system matrix index given electronic state and monomer's state.
                coupled_state_index_in_dimer = search_full_hamiltonian_state_index(coupled_state_exciton_state_index,
                                                                                   coupled_state_monomer1_state_index,
                                                                                   coupled_state_monomer2_state_index);

                if(coupled_state_index_in_dimer == -1){
                    // the state not found. skip to next states.
                    continue;
                }

                franck_condon_factor_monomer1 = d.nonadiabatic_coupled_monomer_state_franck_condon[0][state_monomer1_state_index][j];
                franck_condon_factor_monomer2 = d.nonadiabatic_coupled_monomer_state_franck_condon[1][state_monomer2_state_index][k];

                off_diagonal_matrix_ele = nonadiabatic_coupling * franck_condon_factor_monomer1 * franck_condon_factor_monomer2;

                coupled_state_energy = s.electronic_state_energy[ coupled_state_exciton_state_index ] + monomer1_vib_state_energy_all_pc[coupled_state_monomer1_state_index] + monomer2_vib_state_energy_all_pc[coupled_state_monomer2_state_index];

                state_energy_difference = abs(state_energy - coupled_state_energy);

                if (abs(off_diagonal_matrix_ele) < pow(10,-6) ){
                    ;  // nonadiabatic coupling is too weak, skip that.
                }
                else{
                    if (state_energy_difference == 0){
                        include_nonadiabatic_coupling_bool = true;
                    }
                    else if( abs(off_diagonal_matrix_ele / state_energy_difference) > d.cutoff ){
                        include_nonadiabatic_coupling_bool = true;
                    }
                }

                if (include_nonadiabatic_coupling_bool){
                    // include this off-diagonal coupling.
                    nonadiabatic_off_mat.push_back(off_diagonal_matrix_ele);
                    nonadiabatic_off_irow.push_back(state_index_in_dimer);
                    nonadiabatic_off_icol.push_back(coupled_state_index_in_dimer);
                }

            }
        }

    }

}